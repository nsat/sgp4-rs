use chrono::prelude::*;
use chrono::DateTime;

use thiserror::Error;

mod sgp4_sys;

#[derive(Debug, Error)]
pub enum Error {
    #[error("TLE was malformed: {0}")]
    MalformedTwoLineElement(String),
    #[error("Error in SGP4 propagator")]
    PropagationError,
    #[error("{0}")]
    UnknownError(String),
}

type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
pub struct StateVector {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
}

const TLE_LINE_LENGTH: usize = 69;

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
                "Line 1 is the wrong length. Expected {}, but got {}",
                TLE_LINE_LENGTH,
                line1.len()
            )));
        }

        if line2.len() != TLE_LINE_LENGTH {
            return Err(Error::MalformedTwoLineElement(format!(
                "Line 2 is the wrong length. Expected {}, but got {}",
                TLE_LINE_LENGTH,
                line2.len()
            )));
        }

        let elements = sgp4_sys::to_orbital_elements(
            line1,
            line2,
            sgp4_sys::RunType::Verification,
            sgp4_sys::OperationMode::Improved,
            sgp4_sys::GravitationalConstant::Wgs84,
        )
        .map_err(|e| Error::MalformedTwoLineElement(format!("{:?}", e)))?;

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
        TwoLineElement::new(&lines[0], &lines[1])
    }

    /// Get the epoch of a TwoLineElement.
    pub fn epoch(&self) -> Result<DateTime<Utc>> {
        Ok(self.elements.epoch())
    }

    pub fn propagate_to(&self, t: DateTime<Utc>) -> Result<StateVector> {
        let tle_epoch = self.elements.epoch();
        // TODO: determine correct behaviour for negative prop
        // assert!(t >= tle_epoch);

        let min_since_epoch = (t - tle_epoch).num_minutes() as f64;

        let (r, v) = sgp4_sys::run_sgp4(
            self.elements,
            sgp4_sys::GravitationalConstant::Wgs84,
            min_since_epoch,
        )
        .map_err(|_e| Error::PropagationError)?;

        Ok(StateVector {
            position: r.to_owned(),
            velocity: v.to_owned(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use chrono::Duration;

    #[test]
    fn test_simple_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let s1 = tle.propagate_to(epoch)?;
        let s2 = tle.propagate_to(epoch + Duration::hours(1))?;

        assert_ne!(s1.position, s2.position);
        assert_ne!(s1.velocity, s2.velocity);

        Ok(())
    }

    #[test]
    fn test_negative_time_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let _s1 = tle.propagate_to(epoch);
        let _s2 = tle.propagate_to(epoch - Duration::days(30));

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
}
