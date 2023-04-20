use std::f64::consts::PI;

use chrono::prelude::*;
use chrono::DateTime;
use thiserror::Error;
use uom::si::angle::radian;
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

    pub osculating_elements: OsculatingOrbitalElements,

    pub mean_elements: MeanOrbitalElements,
}

impl StateVector {
    pub fn new(epoch: DateTime<Utc>, position: [f64; 3], velocity: [f64; 3]) -> Self {
        let osculating_elements: OsculatingOrbitalElements = sgp4_sys::to_classical_elements(&position, &velocity).into();
        let mean_elements = osculating_elements.into();
        Self {
            epoch,
            position,
            velocity,
            osculating_elements,
            mean_elements,
        }
    }
}

/// A Keplerian orbital element set.
///
/// This structure contains all of the "classical" orbital elements as derivable from a TLE. We
/// lean on the `uom` crate to provide safe dimensional types which help to avoid bugs related to
/// mixing units of measure.
/// 
/// This elements set is "osculating" in the sense that includes for example the gravitational effect
/// of the Earth's oblate region at the equator. However the osculating elements set cannot be used 
/// to generate a TLE, that requires a mean elements set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OsculatingOrbitalElements {
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

impl From<sgp4_sys::ClassicalOrbitalElements> for OsculatingOrbitalElements {
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

#[cfg(feature = "tlegen")]
impl MeanOrbitalElements {
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

    fn normalize_angle(angle: f64) -> f64 {
        let mut result_angle = angle;
        while result_angle >= 360.0 {
            result_angle -= 360.0;
        }
        while result_angle < 0.0 {
            result_angle += 360.0;
        }
        result_angle
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
        let incl = Self::normalize_angle(self.mean_inclination.get::<angle::degree>());
        let raan = Self::normalize_angle(self.mean_raan.get::<angle::degree>());
        let ecc_int = (self.mean_eccentricity * 10e6).round() as i64;
        let argp = Self::normalize_angle(self.mean_arg_perigee.get::<angle::degree>());
        let ma = Self::normalize_angle(self.mean_mean_anomaly.get::<angle::degree>());
        let consts = sgp4_sys::gravitational_constants();
        let mm = Self::SECONDS_PER_DAY
            / ((2.0 * PI) * (self.mean_semimajor_axis.get::<kilometer>().powi(3) / consts.mu).sqrt());
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeanOrbitalElements {
    mean_semimajor_axis: Length,
    mean_eccentricity: f64,
    mean_inclination: Angle,
    mean_raan: Angle,
    mean_arg_perigee: Angle,
    mean_mean_anomaly: Angle,
}

impl MeanOrbitalElements {
    fn cos(x: f64) -> f64 {
        x.cos()
    }
    fn sin(x: f64) -> f64 {
        x.sin()
    }
    fn tan(x: f64) -> f64 {
        x.tan()
    }
    fn arctan(x: f64) -> f64 {
        x.atan()
    }
    fn arcsin(x: f64) -> f64 {
        x.asin()
    }
    fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }
    fn arctan2(y: f64, x: f64) -> f64 {
        y.atan2(x)
    }
}

// from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html:
// const MAJOR_EARTH_RADIUS: f64 = 6378137.0;  // [m], semi-major axis
// #[allow(non_upper_case_globals)]
// const Rq: f64 = MAJOR_EARTH_RADIUS/1000.0;  // major earth radius in Km
// const J2: f64 = 1.08262668355e-003_f64;


impl From<OsculatingOrbitalElements> for MeanOrbitalElements {
    // suppress these warnings to allow code to be closer to Python original.
    #[allow(non_snake_case, unused_parens)]
    fn from(coe: OsculatingOrbitalElements) -> Self {

        let consts = sgp4_sys::gravitational_constants();

        let Rq: f64 = consts.radiusearthkm;
        let J2: f64 = consts.j2;

        let two_pi = 2.0 * PI;

        let a = coe.semimajor_axis.get::<kilometer>();
        let e = coe.eccentricity;
        let i = coe.inclination.get::<radian>();
        let w = coe.argument_of_perigee.get::<radian>();
        let f = coe.true_anomaly.get::<radian>();
        let u = coe.argument_of_latitude.get::<radian>() % two_pi;
        let OM = coe.raan.get::<radian>();

        let ci = Self::cos(i);
        let E = (2. * Self::arctan(Self::sqrt((1.-e)/(1.+e))*Self::tan(f/2.))) % two_pi;
        let M = E - e * Self::sin(E);
        let gamma2 = -J2/2.*(Rq/a).powi(2);  // the minus means :    Osc ---> Mean
        let eta = Self::sqrt(1.-e.powi(2));
        let gamma2_tag = gamma2/(eta.powi(4));
        let a_r = (1. + e * Self::cos(f)) / (eta.powi(2));
        let a_ave = a + a * gamma2 * ((3.*ci.powi(2)-1.) * (a_r.powi(3) - 1./(eta.powi(3))) + 3.*(1.-ci.powi(2)) * (a_r.powi(3)) * Self::cos(2.*u));
        let de1 = gamma2_tag / 8. * e * (eta.powi(2)) * (1.-11.*ci.powi(2)-40.*(ci.powi(4))/(1.-5.*ci.powi(2))) * Self::cos(2.*w);
        let de = de1 + eta.powi(2) / 2. * (gamma2 * ((3.*ci.powi(2)-1.) / (eta.powi(6)) *
                                           (e*eta + e/(1.+eta) + 3.*Self::cos(f) + 3.*e*(Self::cos(f).powi(2)) +
                                            (e.powi(2))*Self::cos(f).powi(3)) +
                                           3. * (1.-ci.powi(2)) / (eta.powi(6)) * (e + 3.*Self::cos(f) + 3.*e*(Self::cos(f).powi(2)) +
                                                                       (e.powi(2))*Self::cos(f).powi(3)) * Self::cos(2.*u)) -
                                 gamma2_tag * (1.-ci.powi(2)) * (3.*Self::cos(2.*w+f) + Self::cos(2.*w+3.*f)));
        let di = (-e * de1 / ((eta.powi(2))*Self::tan(i)) + gamma2_tag / 2.*ci * Self::sqrt(1.-ci.powi(2)) *
              (3.*Self::cos(2.*w+2.*f) + 3.*e*Self::cos(2.*w+f) + e*Self::cos(2.*w+3.*f)));

        let MWO_ave = (M + w + OM +
                   gamma2_tag / 8. * (eta.powi(3)) * (1.-11.*ci.powi(2)-40.*(ci.powi(4))/(1.-5.*ci.powi(2))) -
                   gamma2_tag / 16. * (2. + e.powi(2)-11.*(2.+3.*e.powi(2))*ci.powi(2) -
                                      40.*(2.+5.*e.powi(2))*(ci.powi(4))/(1.-5.*ci.powi(2)) -
                                      400.*(e.powi(2))*ci.powi(6)/(1.-5.*ci.powi(2)).powi(2)) +
                   gamma2_tag / 4. * (- 6. * (1.-5.*ci.powi(2)) * (f-M+e*Self::sin(f)) +
                                     (3.-5.*ci.powi(2)) * (3.*Self::sin(2.*u)+3.*e*Self::sin(2.*w+f)+e*Self::sin(2.*w+3.*f))) -
                   gamma2_tag / 8. * (e.powi(2)) * ci * (11. + 80.*(ci.powi(2))/(1.-5.*ci.powi(2)) + 200.*(ci.powi(4))/(1.-5.*ci.powi(2)).powi(2)) -
                   gamma2_tag / 2. * ci * (6. * (f-M+e*Self::sin(f)) - 3.*Self::sin(2.*u) - 3.*e*Self::sin(2.*w+f) -
                                          e*Self::sin(2.*w+3.*f)));

        let edM = (gamma2_tag/8.*e*(eta.powi(3))*(1.-11.*ci.powi(2)-40.*(ci.powi(4))/(1.-5.*ci.powi(2))) -
               gamma2_tag/4.*(eta.powi(3))*(2.*(3.*ci.powi(2)-1.)*((a_r*eta).powi(2)+a_r+1.)*Self::sin(f) +
                                      3.*(1.-ci.powi(2))*((-(a_r*eta).powi(2)-a_r+1.)*Self::sin(2.*w+f)+((a_r*eta).powi(2)+a_r+1./3.) *
                                                   Self::sin(2.*w+3.*f))));

        let dOM = (-gamma2_tag/8.*(e.powi(2))*ci*(11.+80.*(ci.powi(2))/(1.-5.*ci.powi(2))+200.*(ci.powi(4))/(1.-5.*ci.powi(2)).powi(2)) -
               gamma2_tag/2.*ci*(6.*(f-M+e*Self::sin(f))-3.*Self::sin(2.*u)-3.*e*Self::sin(2.*w+f)-e*Self::sin(2.*w+3.*f)));

        let d1 = (e+de)*Self::sin(M) + edM*Self::cos(M);
        let d2 = (e+de)*Self::cos(M) - edM*Self::sin(M);
        let M_ave = Self::arctan2(d1, d2) % two_pi;
        let e_ave = Self::sqrt(d1.powi(2) + d2.powi(2));

        let d3 = (Self::sin(i/2.)+Self::cos(i/2.)*di/2.) * Self::sin(OM) + Self::sin(i/2.) * dOM * Self::cos(OM);
        let d4 = (Self::sin(i/2.)+Self::cos(i/2.)*di/2.) * Self::cos(OM) - Self::sin(i/2.) * dOM * Self::sin(OM);
        let OM_ave = Self::arctan2(d3, d4) % two_pi;

        let i_ave = 2. * Self::arcsin(Self::sqrt(d3.powi(2)+d4.powi(2)));

        let w_ave = MWO_ave - M_ave - OM_ave;

        MeanOrbitalElements { 
            mean_semimajor_axis: Length::new::<kilometer>(a_ave),
            mean_eccentricity: e_ave, 
            mean_inclination: Angle::new::<radian>(i_ave), 
            mean_raan: Angle::new::<radian>(OM_ave), 
            mean_arg_perigee: Angle::new::<radian>(w_ave), 
            mean_mean_anomaly: Angle::new::<radian>(M_ave), 
        }
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
    fn test_can_roundtrip_conversion_of_mean_elements_to_tle() -> Result<()> {
        use float_cmp::assert_approx_eq;

        let epoch = Utc.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
        let altitude_km = 408.0;
        let earth_radius_km = 6371.0;
        let coe = MeanOrbitalElements {
            mean_semimajor_axis: Length::new::<kilometer>(altitude_km + earth_radius_km),
            mean_eccentricity: 0.0,
            mean_inclination: Angle::new::<angle::degree>(10.0),
            mean_raan: Angle::new::<angle::degree>(25.0),
            mean_arg_perigee: Angle::new::<angle::degree>(0.0),
            mean_mean_anomaly: Angle::new::<angle::degree>(0.0),
        };
        let tle = coe.as_tle_at(0, epoch);
        println!("Generated TLE:\n{}", tle);
        let sv = TwoLineElement::from_lines(&tle)?.propagate_to(epoch)?;
        let new_coe = sv.mean_elements;
        assert_approx_eq!(
            f64,
            new_coe.mean_semimajor_axis.get::<kilometer>(),
            coe.mean_semimajor_axis.get::<kilometer>(),
            epsilon = 10.0
        );
        assert_approx_eq!(f64, new_coe.mean_eccentricity, coe.mean_eccentricity, epsilon = 0.01);
        Ok(())
    }

    #[test]
    #[cfg(feature = "tlegen")]
    fn test_mean_orbital_elements_to_tle() -> Result<()> {
        use uom::si::length::meter;

        let time = Utc.with_ymd_and_hms(2023, 3, 10, 1, 0, 0).unwrap();

        let coe = MeanOrbitalElements {
            mean_semimajor_axis: Length::new::<meter>(6755925.456114554),
            mean_eccentricity: 0.0013322422991375329,
            mean_inclination: Angle::new::<angle::degree>(0.7850853743058481),
            mean_raan: Angle::new::<angle::degree>(0.4031559142883887),
            mean_arg_perigee: Angle::new::<angle::degree>(2.146362751175218),
            mean_mean_anomaly: Angle::new::<angle::degree>(1.675200832732889),
        };

        let tle_string = coe.as_tle_at(0, time);
        assert_eq!(
            tle_string,
            r#"1 00000U 23001A   23069.04166667  .00000000  00000-0  00000-0 0  9992
2 00000   0.7851   0.4032 0013322   2.1464   1.6752 15.63419485000018"#
        );
        Ok(())
    }
//2 00000  36.9665 -122.4467 0017128 361.6011  -2.1383 14.96542799000016

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
        let mean_elements_1 = svector.mean_elements;
        println!("mean_elements_1:{:?}", mean_elements_1);
        let tle_string = mean_elements_1.as_tle_at(0, epoch);
        let svector_2 = TwoLineElement::from_lines(&tle_string)?.propagate_to(epoch)?;
        let r_2 = svector_2.position;
        let v_2 = svector_2.velocity;
        let mean_elements_2 = svector_2.mean_elements;
        println!("mean_elements_2:{:?}", mean_elements_2);
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
