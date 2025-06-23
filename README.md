# Decoding AGN Feedback with X-arithmetic: From Morphology to Physical Mechanisms

**Authors:**  
Hannah McCall, Irina Zhuravleva, Eugene Churazov

# Usage

This is a code snippet from the paper:
> _Decoding AGN Feedback with X-arithmetic: From Morphology to Physical Mechanisms_.  
> Hannah McCall, Irina Zhuravleva, Eugene Churazov, Congyao Zhang, William Forman, Christine Jones, Yuan Li.  
> ApJ, 2025.

The code applies the X-arithmetic method to produce three map types:
- **`_noshocks`**: excludes weak shocks and sound waves  
- **`_nobubbles`**: excludes bubbles inflated by jets  
- **`_noisobaric`**: excludes slow gas motions and gas cooling  

Users should provide a soft band and a hard band residual image, meaning an image divided by a global model. For best results, the images should have point sources removed. The best-fit spherical or elliptical $\beta$-model parameters should be provided to perform the projection correction. Note that the default emissivity parameters correspond to A2052, and should be calculated and provided for other objects.

## Installation
This code runs in any python environment with Python version 3.11 or later. You can install the requirements via
```
pip install -r requirements.txt
```

## Example

This example uses the galaxy cluster Abell 2052.


```python
soft_resid = fits.getdata("A2052_soft_resid.fits") # Provide the residual images with point sources removed for your cluster
hard_resid = fits.getdata("A2052_hard_resid.fits")

hard_resid = projection_correction(
    hard_resid,
    soft_rc=0.49,  # Example values, replace with best-fit parameters for your cluster
    soft_beta=0.52,
    soft_theta=0.78,
    soft_eps=-0.17,
    hard_rc=0.74,
    hard_beta=0.52,
    hard_theta=0.875,
    hard_eps=-0.2,
    model_center=(122, 122),  # Example center coordinates
    cdelt1=0.000273,  # Example pixel scale in degrees
)

coefficients = calculate_coefficients() # For other clusters, provide the emissivity parameters

xarithmetic_result = calculate_xarithmetic(
    [soft_resid, hard_resid], coefficients
)
```