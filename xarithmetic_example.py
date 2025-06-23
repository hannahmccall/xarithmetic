import numpy as np
from astropy.io import fits
from itertools import product
from scipy.special import gamma

def projection_correction(
    residual: np.ndarray,
    soft_rc: float,
    soft_beta: float,
    soft_theta: float,
    soft_eps: float,
    hard_rc: float,
    hard_beta: float,
    hard_theta: float,
    hard_eps: float,
    model_center: tuple[float, float],
    cdelt1: float,
) -> np.ndarray:
    """Calculates the projection correction for the hard band residual image.
    Parameters
    ----------
    residual: np.ndarray
        *Hard* residual. This will be the residual that will be readjust.
    soft_rc : float
        Beta model (soft band) core radius in arcminutes.
    soft_beta : float
        Beta model (soft band) beta parameter.
    soft_theta : float
        Beta model (soft band) rotation angle in radians.
    soft_eps: float
        Beta model (soft band) ellipticity.
    hard_rc : float
        Beta model (hard band) core radius in arcminutes.
    hard_beta : float
        Beta model (hard band) beta parameter.
    hard_theta : float
        Beta model (hard band) rotation angle in radians.
    hard_eps : float
        Beta model (hard band) ellipticity.
    model_center : tuple[float, float]
        Center in pixel coordinates. Has to be the same for both hard and soft band.
    cdelt1 : float
        Pixel scale in degrees.

    Returns
    -------
    numpy.ndarray, shape=(width, height)
        Corrected residual image.
    """
    mx, my = np.meshgrid(np.arange(residual.shape[0]), np.arange(residual.shape[1]))

    corrs = 1.0 / (soft_rc * gamma(3.0 * soft_beta - 0.5) / gamma(3.0 * soft_beta))
    corrh = 1.0 / (hard_rc * gamma(3.0 * hard_beta - 0.5) / gamma(3.0 * hard_beta))

    mx_eps_soft = (mx - model_center[0]) * np.cos(soft_theta) + (my - model_center[1]) * np.sin(soft_theta)
    my_eps_soft = -(mx - model_center[0]) * np.sin(soft_theta) + (my - model_center[1]) * np.cos(soft_theta)

    mx_eps_hard = (mx - model_center[0]) * np.cos(hard_theta) + (my - model_center[1]) * np.sin(hard_theta)
    my_eps_hard = -(mx - model_center[0]) * np.sin(hard_theta) + (my - model_center[1]) * np.cos(hard_theta)

    rproj_s = abs(
        (
            np.sqrt(mx_eps_soft**2 * (1 - soft_eps) ** 2 + my_eps_soft**2)
            / (1 - soft_eps)
        )
        * cdelt1
        * 60.0
    )
    rproj_h = abs(
        (
            np.sqrt(mx_eps_hard**2 * (1 - hard_eps) ** 2 + my_eps_hard**2)
            / (1 - hard_eps)
        )
        * cdelt1
        * 60.0
    )
    ws = corrs / (1.0 + (rproj_s / soft_rc) ** 2) ** 0.5
    wh = corrh / (1.0 + (rproj_h / hard_rc) ** 2) ** 0.5
    out = residual * ws / wh
    return out


def signed_absmax(x: np.ndarray, axis: int) -> np.ndarray:
    """Returns the element with the largest absolute value along axis of an array."""
    amax = x.max(axis)
    amin = x.min(axis)
    return np.where(-amin > amax, amin, amax)


def calculate_coefficients(
    temperature_params: list[float] = None, emissivity_params: list[float] = None
) -> np.ndarray:
    """Calculates x-arithmetic weight coefficients.

    Parameters
    ----------
    temperature_params : list of floats
        Temperature parameters, one for each perturbation type.
    emissivity_params : list of floats
        Emissivitiy parameters, one for each energy band.

    Returns
    -------
    numpy.ndarray, shape=(n_processes, n_bands)
        Coefficients for x-arithmetic
    """
    temperature_params = temperature_params or [0.66667, 0.0, -1.0]
    emissivity_params = emissivity_params or [-0.5948, 2.3326]

    amp = 2 + np.array(temperature_params)[None] * np.array(emissivity_params)[:, None]
    amp_norm = signed_absmax(
        amp[0][None] * (amp[1] / amp[0])[:, None] - amp[1][None], axis=-1
    )
    return np.stack([(amp[1] / amp[0]) / amp_norm, -1.0 / amp_norm], axis=-1)


def calculate_xarithmetic(
    residuals: list[np.ndarray], coefficients: np.ndarray
) -> np.ndarray:
    """Main routine.

    Parameters
    ----------
    residuals : list[numpy.ndarray], shape=(width, height)
        Residuals of the images for each band.
    coefficients : numpy.ndarray, shape=(n_processes, n_bands)
        X-arthmetic coefficients.

    Returns
    -------
    numpy.ndarray, shape=(n_processes, width, height)
        Processed images.
    """
    return np.sum(np.array(residuals)[None] * coefficients[..., None, None], axis=1)