//! Fixture for mantra-testrun's self-check. See reqs/req.md for what each
//! requirement is expected to resolve to and why.
pub mod elsewhere;
pub mod never;

#[mantra_macros::req_satisfied("p-1")]
pub fn one() -> u32 {
    1
}

#[mantra_macros::req_satisfied("p-3")]
pub fn three() -> u32 {
    3
}

#[mantra_macros::req_satisfied("p-5")]
pub fn five() -> u32 {
    5
}

#[mantra_macros::req_satisfied("p-7", "p-8")]
pub fn seven_eight() -> u32 {
    78
}

#[mantra_macros::req_satisfied("p-9")]
pub fn nine() -> u32 {
    9
}

/// Trait declarations have no executable lines: a satisfy trace here can never
/// be confirmed by coverage (p-11).
#[mantra_macros::req_satisfied("p-11")]
pub trait Eleven {
    fn eleven(&self) -> u32;
}

pub struct E;

impl Eleven for E {
    fn eleven(&self) -> u32 {
        11
    }
}

#[mantra_macros::req_satisfied("p-12")]
pub fn twelve() -> u32 {
    12
}

#[cfg(test)]
mod tests {
    use super::*;

    #[mantra_macros::req_test("p-1")]
    #[test]
    fn p1_executes_one() {
        assert_eq!(one(), 1);
    }

    #[mantra_macros::req_test("p-2")]
    #[test]
    fn p2_passes_without_executing_never_executed() {
        assert_eq!(1 + 1, 2);
    }

    #[mantra_macros::req_test("p-3")]
    #[test]
    fn p3_fails_on_purpose() {
        assert_eq!(three(), 4, "deliberate failure for the self-check");
    }

    #[mantra_macros::req_test("p-6")]
    #[test]
    fn p6_traced_but_does_not_execute_six() {
        assert_eq!(elsewhere::helper(), 0);
    }

    #[test]
    fn p6_untraced_executes_six() {
        assert_eq!(elsewhere::six(), 6);
    }

    #[mantra_macros::req_test("p-7")]
    #[test]
    fn p7_first_verifier() {
        assert_eq!(seven_eight(), 78);
    }

    #[mantra_macros::req_test("p-7", "p-8")]
    #[test]
    fn p7_and_p8_second_verifier() {
        assert_eq!(seven_eight(), 78);
    }

    #[mantra_macros::req_test("p-9")]
    #[test]
    #[ignore = "skipped on purpose for the self-check"]
    fn p9_ignored() {
        assert_eq!(nine(), 9);
    }

    #[mantra_macros::req_test("p-11")]
    #[test]
    fn p11_calls_trait_method() {
        assert_eq!(E.eleven(), 11);
    }

    #[mantra_macros::req_test("p-12")]
    #[test]
    fn p12_executes_twelve() {
        assert_eq!(twelve(), 12);
    }

    #[mantra_macros::req_test("p-12")]
    #[test]
    fn p12_does_not_execute_twelve() {
        assert_eq!(2 + 2, 4);
    }
}
