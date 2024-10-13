use chrono::prelude::*;
use chrono::DateTime;
use thiserror::Error;
use uom::si::{angle, angular_velocity::radian_per_second, f64::*, length::kilometer};

mod sgp4_sys;
#[cfg(feature = "tlegen")]
mod tlegen;

#[derive(Debug, Error, PartialEq)]
pub enum Error {
    #[error("TLE was malformed: {0}")]
    MalformedTwoLineElement(String),
    #[error("{0}")]
    UnknownError(String),
    #[error(transparent)]
    PropagationError(#[from] sgp4_sys::Error),
    #[error("Optimization error: {0}")]
    OptimizationError(String),
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
    fn test_epoch_has_subsecond_precision() {
        let lines = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
                     2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";
        let tle = TwoLineElement::from_lines(lines).unwrap();
        assert_eq!(tle.epoch().unwrap().nanosecond(), 452_818_000);
        let t = NaiveDate::from_ymd_opt(2020, 5, 29)
            .unwrap()
            .and_hms_micro_opt(1, 2, 3, 0)
            .unwrap()
            .and_utc();
        let position = tle.propagate_to(t).unwrap().position;
        assert_eq!((position[0] * 1e3) as i64, 4_262_298);
        assert_eq!((position[1] * 1e3) as i64, 403_136);
        assert_eq!((position[2] * 1e3) as i64, -5_284_434);
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
}
