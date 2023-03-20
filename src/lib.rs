use chrono::prelude::*;
use chrono::DateTime;
use thiserror::Error;
use uom::si::{angle, angular_velocity::radian_per_second, f64::*, length::kilometer};

mod sgp4_sys;

#[derive(Debug, Error, PartialEq)]
pub enum Error {
    #[error("TLE was malformed: {0}")]
    MalformedTwoLineElement(String),
    #[error("{0}")]
    UnknownError(String),
    #[error(transparent)]
    PropagationError(#[from] sgp4_sys::Error),
}

type Result<T> = std::result::Result<T, Error>;

/// A state vector, in the TEME-ECI coordinate frame, for an orbiting body.
///
/// To obtain the state of an object at a specific time, use the propagation functions provided by
/// [TwoLineElement].
#[derive(Debug, Clone, Copy)]
pub struct StateVector {
    pub epoch: DateTime<Utc>,

    /// The satellite position in km.
    pub position: [f64; 3],

    /// The satellite velocity in km/s.
    pub velocity: [f64; 3],

    pub coe: ClassicalOrbitalElements,
}

impl StateVector {
    pub fn new(epoch: DateTime<Utc>, position: [f64; 3], velocity: [f64; 3]) -> Self {
        Self {
            epoch,
            position,
            velocity,
            coe: sgp4_sys::to_classical_elements(&position, &velocity).into(),
        }
    }

    pub fn semilatus_rectum(&self) -> Length {
        self.coe.semilatus_rectum
    }

    pub fn semimajor_axis(&self) -> Length {
        self.coe.semimajor_axis
    }

    pub fn inclination(&self) -> Angle {
        self.coe.inclination
    }

    pub fn raan(&self) -> Angle {
        self.coe.raan
    }

    pub fn mean_anomaly(&self) -> Angle {
        self.coe.mean_anomaly
    }

    pub fn true_anomaly(&self) -> Angle {
        self.coe.true_anomaly
    }

    pub fn eccentricity(&self) -> f64 {
        self.coe.eccentricity
    }

    pub fn longitude_of_periapsis(&self) -> Angle {
        self.coe.longitude_of_periapsis
    }

    pub fn true_longitude(&self) -> Angle {
        self.coe.true_longitude
    }

    pub fn argument_of_perigee(&self) -> Angle {
        self.coe.argument_of_perigee
    }

    pub fn argument_of_latitude(&self) -> Angle {
        self.coe.argument_of_latitude
    }
}

/// A Keplerian orbital element set.
///
/// This structure contains all of the "classical" orbital elements as derivable from a TLE. We
/// lean on the `uom` crate to provide safe dimensional types which help to avoid bugs related to
/// mixing units of measure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassicalOrbitalElements {
    pub semilatus_rectum: Length,
    pub semimajor_axis: Length,
    pub eccentricity: f64,
    pub inclination: Angle,
    pub raan: Angle,
    pub argument_of_perigee: Angle,
    pub true_anomaly: Angle,
    pub mean_anomaly: Angle,
    pub argument_of_latitude: Angle,
    pub true_longitude: Angle,
    pub longitude_of_periapsis: Angle,
}

impl From<sgp4_sys::ClassicalOrbitalElements> for ClassicalOrbitalElements {
    fn from(coe: sgp4_sys::ClassicalOrbitalElements) -> Self {
        let semilatus_rectum = Length::new::<kilometer>(coe.p);
        let semimajor_axis = Length::new::<kilometer>(coe.a);
        let inclination = Angle::new::<angle::radian>(coe.incl);
        let raan = Angle::new::<angle::radian>(coe.omega);
        let mean_anomaly = Angle::new::<angle::radian>(coe.m);
        let true_anomaly = Angle::new::<angle::radian>(coe.nu);
        let eccentricity = coe.ecc;
        let longitude_of_periapsis = Angle::new::<angle::radian>(coe.lonper);
        let true_longitude = Angle::new::<angle::radian>(coe.truelon);
        let argument_of_perigee = Angle::new::<angle::radian>(coe.argp);
        let argument_of_latitude = Angle::new::<angle::radian>(coe.arglat);

        Self {
            semilatus_rectum,
            semimajor_axis,
            eccentricity,
            inclination,
            raan,
            argument_of_perigee,
            true_anomaly,
            mean_anomaly,
            argument_of_latitude,
            true_longitude,
            longitude_of_periapsis,
        }
    }
}

impl From<StateVector> for ClassicalOrbitalElements {
    fn from(sv: StateVector) -> Self {
        sv.coe
    }
}

#[cfg(feature = "tlegen")]
impl ClassicalOrbitalElements {
    const SECONDS_PER_DAY: f64 = 24.0 * 60.0 * 60.0;

    /// Create a formatted Two Line Element string from a Keplerian orbital element set for testing
    /// purposes.
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
    pub fn as_tle_at(&self, catalog_num: u8, epoch: DateTime<Utc>) -> String {
        format!(
            "{}\n{}",
            self.tle_line_1(catalog_num, epoch),
            self.tle_line_2(catalog_num)
        )
    }

    #[cfg(feature = "tlegen")]
    fn tle_line_1(&self, catalog_num: u8, epoch: DateTime<Utc>) -> String {
        let epoch_year = epoch.year() % 100;
        let epoch_day = epoch.ordinal();
        let epoch_day_fraction = epoch.num_seconds_from_midnight() as f64 / Self::SECONDS_PER_DAY;
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
        Self::add_tle_checksum(line)
    }

    fn tle_line_2(&self, catalog_num: u8) -> String {
        use std::f64::consts::PI;

        let incl = self.inclination.get::<angle::degree>();
        let raan = self.raan.get::<angle::degree>();
        let ecc_int = (self.eccentricity * 10e6).round() as i64;
        let argp = self.argument_of_perigee.get::<angle::degree>();
        let ma = self.mean_anomaly.get::<angle::degree>();
        let consts = sgp4_sys::gravitational_constants();
        let mm = Self::SECONDS_PER_DAY
            / ((2.0 * PI) * (self.semimajor_axis.get::<kilometer>().powi(3) / consts.mu).sqrt());
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
        Self::add_tle_checksum(line)
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
}

const TLE_LINE_LENGTH: usize = 69;

/// A parsed, valid Two Line Element data set which can be used for orbital propagation.
///
/// Internally this uses SGP4's own structure representation. Various fields which are useful for
/// analysis or simulation are exposed via methods like `raan()`/`set_raan()`
/// which allow access to the values and direct modification of the underlying orbital element set
/// in a type-safe manner. The `uom` crate provides dimensional analysis to help avoid
/// unit-of-measure errors which can otherwise be quite difficult to detect.
#[derive(Clone)]
pub struct TwoLineElement {
    elements: sgp4_sys::OrbitalElementSet,
}

impl TwoLineElement {
    /// Create a validated TwoLineElement from a string.
    pub fn new(line1: &str, line2: &str) -> Result<TwoLineElement> {
        let line1 = line1.trim();
        let line2 = line2.trim();

        if line1.len() != TLE_LINE_LENGTH {
            return Err(Error::MalformedTwoLineElement(format!(
                "Line 1 is the wrong length. Expected {}, but got {}\n{}",
                TLE_LINE_LENGTH,
                line1.len(),
                line1
            )));
        }

        if line2.len() != TLE_LINE_LENGTH {
            return Err(Error::MalformedTwoLineElement(format!(
                "Line 2 is the wrong length. Expected {}, but got {}\n{}",
                TLE_LINE_LENGTH,
                line2.len(),
                line2
            )));
        }

        let elements = sgp4_sys::to_orbital_elements(
            line1,
            line2,
            sgp4_sys::RunType::Verification,
            sgp4_sys::OperationMode::Improved,
            sgp4_sys::GravitationalConstant::Wgs84,
        )
        .map_err(|e| Error::MalformedTwoLineElement(e.to_string()))?;

        Ok(TwoLineElement { elements })
    }

    /// Create a TwoLineElement from a string containing both lines, and optionally a header line.
    pub fn from_lines(combined_lines: &str) -> Result<TwoLineElement> {
        let lines: Vec<_> = {
            let mut ls: Vec<_> = combined_lines
                .split('\n')
                .filter(|s| !s.is_empty())
                .collect();
            if ls.len() == 3 {
                ls.split_off(1)
            } else if ls.len() == 2 {
                ls
            } else {
                return Err(Error::MalformedTwoLineElement(format!(
                    "Expected two lines, got {}",
                    ls.len()
                )));
            }
        };
        TwoLineElement::new(lines[0], lines[1])
    }

    /// Get the epoch of a TwoLineElement.
    pub fn epoch(&self) -> Result<DateTime<Utc>> {
        Ok(self.elements.epoch())
    }

    pub fn mean_motion(&self) -> AngularVelocity {
        AngularVelocity::new::<radian_per_second>(self.elements.mean_motion() / 60.)
    }

    /// Propagate a TwoLineElement to the given time to obtain a state vector for the object.
    pub fn propagate_to(&self, t: DateTime<Utc>) -> Result<StateVector> {
        let tle_epoch = self.elements.epoch();
        // TODO: determine correct behaviour for negative prop
        // assert!(t >= tle_epoch);

        let min_since_epoch = (t - tle_epoch).num_milliseconds() as f64 / 60_000.;

        let (r, v) = sgp4_sys::run_sgp4(
            self.elements,
            sgp4_sys::GravitationalConstant::Wgs84,
            min_since_epoch,
        )?;

        Ok(StateVector::new(t, r.to_owned(), v.to_owned()))
    }
}

/// Wrapper type representing a Julian day.
///
/// This is the number of days since the start of the Julian astronomical calendar in 4713 BC, used
/// to provide a consistent time reference for astronomical calculations.
pub struct JulianDay(f64);

impl From<DateTime<Utc>> for JulianDay {
    fn from(d: DateTime<Utc>) -> Self {
        JulianDay(sgp4_sys::datetime_to_julian_day(d))
    }
}

impl From<JulianDay> for DateTime<Utc> {
    fn from(jd: JulianDay) -> Self {
        sgp4_sys::julian_day_to_datetime(jd.0)
    }
}

/// Wrapper type representing the angular form of Greenwich Mean Sidereal Time.
///
/// This is primarily used to account for the Earth's rotation during conversion between fixed and
/// inertial coordinate frames. Note that this is an angle measured in radians, and not a "time" as
/// such. The value may range from 0 to 2π.
pub struct GreenwichMeanSiderealTime(f64);

impl GreenwichMeanSiderealTime {
    pub fn as_radians(&self) -> f64 {
        self.0
    }
}

impl From<DateTime<Utc>> for GreenwichMeanSiderealTime {
    fn from(d: DateTime<Utc>) -> Self {
        GreenwichMeanSiderealTime(sgp4_sys::datetime_to_gstime(d))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use chrono::Duration;
    use float_cmp::approx_eq;

    fn vecs_eq(l: &[f64; 3], r: &[f64; 3]) -> bool {
        approx_eq!(f64, l[0], r[0]) && approx_eq!(f64, l[1], r[1]) && approx_eq!(f64, l[2], r[2])
    }

    #[test]
    fn test_simple_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let s1 = tle.propagate_to(epoch)?;
        let s2 = tle.propagate_to(epoch + Duration::hours(1))?;

        assert!(!vecs_eq(&s1.position, &s2.position));
        assert!(!vecs_eq(&s1.velocity, &s2.velocity));

        Ok(())
    }

    #[test]
    fn test_decay_error() -> Result<()> {
        let line1 = "1 43051U 17071Q   22046.92182028  .07161566  12340-4  74927-3 0  9993";
        let line2 = "2 43051  51.6207 236.5853 0009084 284.2762  75.7254 16.36736354237455";
        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let s2 = tle.propagate_to(epoch + Duration::days(30));
        assert!(s2.is_err());
        assert_eq!(
            s2.unwrap_err(),
            Error::PropagationError(sgp4_sys::Error::SatelliteDecay)
        );
        Ok(())
    }

    #[test]
    fn mean_motion_round_trip() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";
        let tle = TwoLineElement::new(line1, line2)?;
        let mean_motion = tle.mean_motion().get::<radian_per_second>() * 24. * 60.0 * 60.0
            / (2. * std::f64::consts::PI);
        assert!(approx_eq!(f64, mean_motion, 15.493968, epsilon = 0.01));
        Ok(())
    }

    #[test]
    fn test_negative_time_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let s1 = tle.propagate_to(epoch)?;
        let s2 = tle.propagate_to(epoch - Duration::days(30))?;

        assert!(!vecs_eq(&s1.position, &s2.position));
        assert!(!vecs_eq(&s1.velocity, &s2.velocity));

        Ok(())
    }

    #[test]
    fn test_tle_from_lines() -> Result<()> {
        let lines = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
                     2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let _tle = TwoLineElement::from_lines(lines)?;
        Ok(())
    }

    #[test]
    fn test_tle_from_lines_with_header() -> Result<()> {
        let lines = "ISS (ZARYA)
                     1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
                     2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let _tle = TwoLineElement::from_lines(lines)?;
        Ok(())
    }

    #[test]
    fn test_tle_from_lines_with_surrounding_whitespace() -> Result<()> {
        let lines = "\nISS (ZARYA)
                     1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
                     2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767\n";

        let _tle = TwoLineElement::from_lines(lines)?;

        Ok(())
    }

    #[test]
    fn test_julian_day_identity() {
        let t = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        assert_eq!(DateTime::<Utc>::from(JulianDay::from(t)), t);
    }

    #[test]
    fn test_gmst_conversion() {
        let t = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        let a: f64 = 100.1218209532; // GMST for 2020-01-01T00:00:00 in degrees
        let a_rad = a.to_radians();
        assert!(sgp4_sys::close(
            GreenwichMeanSiderealTime::from(t).as_radians(),
            a_rad
        ));
    }

    #[test]
    #[cfg(feature = "tlegen")]
    fn test_can_roundtrip_conversion_of_classical_elements_to_tle() -> Result<()> {
        use float_cmp::assert_approx_eq;

        let epoch = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        let altitude_km = 408.0;
        let earth_radius_km = 6371.0;
        let coe = ClassicalOrbitalElements {
            semilatus_rectum: Length::new::<kilometer>(altitude_km + earth_radius_km),
            semimajor_axis: Length::new::<kilometer>(altitude_km + earth_radius_km),
            eccentricity: 0.0,
            inclination: Angle::new::<angle::degree>(10.0),
            raan: Angle::new::<angle::degree>(25.0),
            argument_of_perigee: Angle::new::<angle::degree>(0.0),
            true_anomaly: Angle::new::<angle::degree>(0.0),
            mean_anomaly: Angle::new::<angle::degree>(0.0),
            argument_of_latitude: Angle::new::<angle::degree>(0.0),
            true_longitude: Angle::new::<angle::degree>(0.0),
            longitude_of_periapsis: Angle::new::<angle::degree>(0.0),
        };
        let tle = coe.as_tle_at(0, epoch);
        println!("Generated TLE:\n{}", tle);
        let sv = TwoLineElement::from_lines(&tle)?.propagate_to(epoch)?;
        let new_coe = ClassicalOrbitalElements::from(sv);
        assert_approx_eq!(
            f64,
            new_coe.semimajor_axis.get::<kilometer>(),
            coe.semimajor_axis.get::<kilometer>(),
            epsilon = 10.0
        );
        assert_approx_eq!(f64, new_coe.eccentricity, coe.eccentricity, epsilon = 0.01);
        Ok(())
    }

    #[test]
    #[cfg(feature = "tlegen")]
    fn test_classical_orbital_elements_to_tle() -> Result<()> {
        use uom::si::length::meter;

        let time = Utc.with_ymd_and_hms(2023, 3, 10, 1, 0, 0).unwrap();

        let coe = ClassicalOrbitalElements {
            semilatus_rectum: Length::new::<meter>(6755913.465228223),
            semimajor_axis: Length::new::<meter>(6755925.456114554),
            eccentricity: 0.0013322422991375329,
            inclination: Angle::new::<angle::degree>(0.7850853743058481),
            raan: Angle::new::<angle::degree>(0.4031559142883887),
            argument_of_perigee: Angle::new::<angle::degree>(2.146362751175218),
            true_anomaly: Angle::new::<angle::degree>(1.6778503457504903),
            mean_anomaly: Angle::new::<angle::degree>(1.675200832732889),
            argument_of_latitude: Angle::new::<angle::degree>(999999.1),
            true_longitude: Angle::new::<angle::degree>(999999.1),
            longitude_of_periapsis: Angle::new::<angle::degree>(999999.1),
        };

        let tle_string = coe.as_tle_at(0, time);
        assert_eq!(
            tle_string,
            r#"1 00000U 23001A   23069.04166667  .00000000  00000-0  00000-0 0  9992
2 00000   0.7851   0.4032 0013322   2.1464   1.6752 15.63419485000018"#
        );
        Ok(())
    }

    #[test]
    #[cfg(feature = "tlegen")]
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
        let tle_string = svector.coe.as_tle_at(0, epoch);
        let svector_2 = TwoLineElement::from_lines(&tle_string)?.propagate_to(epoch)?;
        let r_2 = svector_2.position;
        let v_2 = svector_2.velocity;
        println!("r_1:{:?}", r_1);
        println!("r_2:{:?}", r_2);
        println!("v_1:{:?}", v_1);
        println!("v_2:{:?}", v_2);
        assert_approx_eq!(f64, r_1[0], r_2[0], epsilon = 10.);
        assert_approx_eq!(f64, r_1[1], r_2[1], epsilon = 10.);
        assert_approx_eq!(f64, r_1[2], r_2[2], epsilon = 10.);
        assert_approx_eq!(f64, v_1[0], v_2[0], epsilon = 0.01);
        assert_approx_eq!(f64, v_1[1], v_2[1], epsilon = 0.01);
        assert_approx_eq!(f64, v_1[2], v_2[2], epsilon = 0.01);
        Ok(())
    }
}
