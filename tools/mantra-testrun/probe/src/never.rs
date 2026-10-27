/// Satisfies p-2 but is never executed by any test.
#[mantra_macros::req_satisfied("p-2")]
pub fn never_executed() -> u32 {
    2
}
