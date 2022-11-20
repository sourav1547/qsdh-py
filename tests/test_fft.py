from pypairing import ZR, G1, blsfft


def test_fft():
    # Given
    coeffs = [0, 1]
    p = 13
    omega = 5
    n = 4

    # When
    fft_rep = blsfft(coeffs, omega, n)

    # Then
    assert fft_rep == [1, 5, 12, 8]
    