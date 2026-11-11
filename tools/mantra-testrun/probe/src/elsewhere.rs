/// Satisfies p-6; executed only by a test that carries no `req_test`.
#[mantra_macros::req_satisfied("p-6")]
pub fn six() -> u32 {
    6
}

/// Executed by the test that claims to verify p-6, but does not satisfy it.
pub fn helper() -> u32 {
    0
}
