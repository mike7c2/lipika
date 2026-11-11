#[mantra_macros::req_test("p-5")]
#[test]
fn p5_verified_from_integration_test() {
    assert_eq!(probe::five(), 5);
}
