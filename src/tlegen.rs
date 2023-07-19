use argmin::{
    core::{CostFunction, Executor, observers::{SlogLogger, ObserverMode}},
    solver::neldermead::NelderMead,
};
use chrono::{DateTime, Datelike, Timelike, Utc};
use uom::{
    si::{
        angle::{self, degree},
        f64::{Angle, Length},
        length::kilometer,
    },
    ConstZero,
};

use crate::{sgp4_sys, ClassicalOrbitalElements, Error, Result, StateVector, TwoLineElement};

const SECONDS_PER_DAY: f64 = 24.0 * 60.0 * 60.0;

impl ClassicalOrbitalElements {
    pub fn as_tle_at(&self, catalog_num: u8, epoch: DateTime<Utc>) -> String {
        let tle = format!(
            "{}\n{}",
            tle_line_1(catalog_num, epoch),
            tle_line_2(
                catalog_num,
                self.inclination,
                self.raan,
                self.eccentricity,
                self.argument_of_perigee,
                self.mean_anomaly,
                self.semimajor_axis
            )
        );
        tle
    }
}

impl StateVector {
    /// Find a TLE string that propagates to the state vector at a given epoch
    ///
    /// Note that the generated TLE has the following simplifications:
    /// 1. It assumes that the epoch and the launch date are the same.
    /// 2. Launch number is assumed to be 1, and the launch piece is A.
    /// 3. Element set number is always 999.
    /// 4. Mean motion derivatives and ballistic coefficient are set to zero.
    /// 5. The orbit number is assumed to be zero.
    ///
    /// Because of these simplifications, the elements of the generated TLE are not guaranteed to
    /// exactly match those of the original element set. This function should not be used for
    /// production applications.
    pub fn as_tle_at(&self, catalog_num: u8, epoch: DateTime<Utc>) -> Result<String> {
        // The orbital elements associated with the state vector are osculating/instantaneous
        // whereas the TLE must be based on mean elements. To find a TLE which propagates to the
        // required cartesian state vector we use a numerical optimization approach, based on the
        // Nelder-Mead algorithm. We use the six parameters of the instantaneous orbital elements to
        // define the initial simplex. For the objective function we propagate the TLE and calculate
        // the error against the target position / velocity using sum-of-squares.
        let cost = FindTleProblem {
            epoch,
            position: self.position,
            velocity: self.velocity,
        };
        let init_param: Vec<f64> = vec![
            self.coe.inclination.get::<degree>(),
            self.coe.raan.get::<degree>(),
            self.coe.eccentricity,
            self.coe.argument_of_perigee.get::<degree>(),
            self.coe.mean_anomaly.get::<degree>(),
            self.coe.semimajor_axis.get::<kilometer>(),
        ];
        let mut initial_simplex: Vec<Vec<f64>> = vec![];
        for _i in 0..7 {
            // initial simplex requires n+1 vertices
            initial_simplex.push(init_param.clone());
        }
        let perturbations = [0.1, 0.1, 0.01, 1.0, 1.0, 1.0];
        for i in 0..6 {
            // use a custom offset for each parameter
            initial_simplex[i][i] += perturbations[i];
        }
        let solver: NelderMead<Vec<f64>, f64> = NelderMead::new(initial_simplex)
        .with_alpha(0.9).expect("error")
        .with_gamma(1.1).expect("error")
        .with_rho(0.25).expect("error")
        .with_sigma(0.1).expect("error");
        let res = Executor::new(cost, solver)
            .configure(|state| state.param(init_param).max_iters(1000).target_cost(0.0))
            //.add_observer(SlogLogger::term(), ObserverMode::Always)
            .run();
        match res {
            Ok(opt_res) => {
                let best_param = opt_res.state().best_param.as_ref().unwrap();
                let tle = format!(
                    "{}\n{}",
                    tle_line_1(catalog_num, epoch),
                    params_to_tle_line2(catalog_num, best_param)
                );
                Ok(tle)
            }
            Err(opt_err) => Err(Error::OptimizationError(opt_err.to_string())),
        }
    }
}

fn tle_line_1(catalog_num: u8, epoch: DateTime<Utc>) -> String {
    let epoch_year = epoch.year() % 100;
    let epoch_day = epoch.ordinal();
    let epoch_day_fraction = epoch.num_seconds_from_midnight() as f64 / SECONDS_PER_DAY;
    let epoch_day_fraction_int = (epoch_day_fraction * 100000000.0).round() as i64;
    let line = format!(
        "1 {0:05}U {1:2}001A   {1:2}{2:03}.{3:08}  .00000000  00000-0  00000-0 0  999",
        // |-----| |---------| |---||---| |-----| |--------| |------| |------| ^ |--|
        // 3-8     10-17       19      23 25-32   34-43      45-52    54-61      65 68
        catalog_num,
        epoch_year,
        epoch_day,
        epoch_day_fraction_int
    );
    add_tle_checksum(line)
}

fn tle_line_2(
    catalog_num: u8,
    inclination: Angle,
    raan: Angle,
    eccentricity: f64,
    argument_of_perigee: Angle,
    mean_anomaly: Angle,
    semimajor_axis: Length,
) -> String {
    use std::f64::consts::PI;

    let incl = normalize_angle(inclination).get::<angle::degree>();
    let raan = normalize_angle(raan).get::<angle::degree>();
    let ecc_int = (eccentricity * 10e6).round() as i64;
    let argp = normalize_angle(argument_of_perigee).get::<angle::degree>();
    let ma = normalize_angle(mean_anomaly).get::<angle::degree>();
    let consts = sgp4_sys::gravitational_constants();
    let mm = SECONDS_PER_DAY
        / ((2.0 * PI) * (semimajor_axis.get::<kilometer>().powi(3) / consts.mu).sqrt());
    let line = format!(
        "2 {0:05} {1:>8.4} {2:>8.4} {3:07} {4:>8.4} {5:>8.4} {6:>11.8}00001",
        // |----| |------| |------| |----| |------| |------| |-------||---|
        // 3-7    9-16     18-25    27-33  35-42    44-51    53-63    64-68
        catalog_num,
        incl,
        raan,
        ecc_int,
        argp,
        ma,
        mm
    );
    add_tle_checksum(line)
}

fn add_tle_checksum(mut line: String) -> String {
    let checksum = line.chars().fold(0, |acc, c| {
        acc + match c {
            '-' => 1,
            c if c.is_ascii_digit() => c.to_digit(10).unwrap(),
            _ => 0,
        }
    }) % 10;
    line.push_str(&checksum.to_string());
    line
}

struct FindTleProblem {
    pub epoch: DateTime<Utc>,
    pub position: [f64; 3],
    pub velocity: [f64; 3],
}

fn params_to_tle_line2(catalog_num: u8, param: &[f64]) -> String {
    let inclination = Angle::new::<degree>(param[0]);
    let raan = Angle::new::<degree>(param[1]);
    let eccentricity = param[2];
    let argument_of_perigee = Angle::new::<degree>(param[3]);
    let mean_anomaly = Angle::new::<degree>(param[4]);
    let semimajor_axis = Length::new::<kilometer>(param[5]);

    tle_line_2(
        catalog_num,
        normalize_angle(inclination),
        normalize_angle(raan),
        clamp_eccentricity(eccentricity),
        normalize_angle(argument_of_perigee),
        normalize_angle(mean_anomaly),
        semimajor_axis.max(Length::ZERO),
    )
}

fn clamp_eccentricity(ecc: f64) -> f64 {
    ecc.max(0.0).min(1.0)
}

fn normalize_angle(angle: Angle) -> Angle {
    let mut normalized = angle;
    while normalized < Angle::ZERO {
        normalized += Angle::FULL_TURN;
    }
    while normalized >= Angle::FULL_TURN {
        normalized -= Angle::FULL_TURN;
    }
    normalized
}

impl CostFunction for FindTleProblem {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> std::result::Result<f64, argmin::core::Error> {
        let catalog_num = 1;

        let tle_line_1 = tle_line_1(catalog_num, self.epoch);
        let tle_line_2 = params_to_tle_line2(catalog_num, param);
        let tle = TwoLineElement::new(&tle_line_1, &tle_line_2)?;
        let prop_sv = tle.propagate_to(self.epoch)?;

        let error = (self.position[0] - prop_sv.position[0]).powi(2)
            + (self.position[1] - prop_sv.position[1]).powi(2)
            + (self.position[2] - prop_sv.position[2]).powi(2)
            + (self.velocity[0] - prop_sv.velocity[0]).powi(2)
            + (self.velocity[1] - prop_sv.velocity[1]).powi(2)
            + (self.velocity[2] - prop_sv.velocity[2]).powi(2);

        Ok(error)
    }
}

#[cfg(test)]
mod tests {
    use crate::TwoLineElement;

    use super::*;

    use chrono::TimeZone;

    #[test]
    fn test_can_roundtrip_state_vector_plus_epoch_to_tle() -> Result<()> {
        use float_cmp::assert_approx_eq;

        let epoch = Utc.with_ymd_and_hms(2021, 5, 25, 0, 0, 0).unwrap();
        let r_1 = [
            -3767.0783048821595,
            -5832.3746513067335,
            0.013350841794354097,
        ];
        let v_1 = [5.087843659697572, -3.2858873951805836, 4.561428718239809];
        let svector = StateVector::new(epoch, r_1, v_1);
        let tle_string = svector.as_tle_at(0, epoch).unwrap();
        println!("tle_string:\n{}", tle_string);
        let svector_2 = TwoLineElement::from_lines(&tle_string)?.propagate_to(epoch)?;
        let r_2 = svector_2.position;
        let v_2 = svector_2.velocity;
        println!("r_1:{:?}", r_1);
        println!("r_2:{:?}", r_2);
        println!("v_1:{:?}", v_1);
        println!("v_2:{:?}", v_2);
        assert_approx_eq!(f64, r_1[0], r_2[0], epsilon = 0.01);
        assert_approx_eq!(f64, r_1[1], r_2[1], epsilon = 0.01);
        assert_approx_eq!(f64, r_1[2], r_2[2], epsilon = 0.01);
        assert_approx_eq!(f64, v_1[0], v_2[0], epsilon = 0.01);
        assert_approx_eq!(f64, v_1[1], v_2[1], epsilon = 0.01);
        assert_approx_eq!(f64, v_1[2], v_2[2], epsilon = 0.01);
        Ok(())
    }


    #[test]
    fn test_roundtrip_tle_to_tle() -> Result<()> {
        let epoch = Utc.with_ymd_and_hms(2021, 5, 25, 0, 0, 0).unwrap();
        let tle_1 = "1 00000U 21001A   21145.00000000  .00000000  00000-0  00000-0 0  9997\n2 00000  36.9006 237.1418 0013279   1.4043 318.6732 14.97334669000013".to_string();
        println!("tle_1:\n{}", tle_1);
        let svector = TwoLineElement::from_lines(&tle_1)?.propagate_to(epoch)?;
        let tle_2 = svector.as_tle_at(0, epoch).unwrap();
        println!("tle_2:\n{}", tle_2);
        Ok(())
    }

    #[test]
    fn test_roundtrip_tle_to_tle_2() -> Result<()> {
        let epoch = Utc.with_ymd_and_hms(2021, 5, 25, 0, 0, 0).unwrap();
        let tle_1 = "1 00000U 21001A   21145.00000000  .00000000  00000-0  00000-0 0  9997\n2 00000  36.9144 237.1225 1121181   3.5239 316.6354 14.96118753000010".to_string();
        println!("tle_1:\n{}", tle_1);
        let svector = TwoLineElement::from_lines(&tle_1)?.propagate_to(epoch)?;
        let tle_2 = svector.as_tle_at(0, epoch).unwrap();
        println!("tle_2:\n{}", tle_2);
        Ok(())
    }
}
