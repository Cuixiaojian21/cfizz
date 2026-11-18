# cfizz branding assets

This directory contains the cfizz logo and the GitHub social preview image.

## Files

| File | Use |
|---|---|
| `cfizz-logo.png` | original 1024 × 1024 logo on white |
| `cfizz-logo-transparent.png` | original 1024 × 1024 logo with transparent background |
| `cfizz-logo-256.png` | 256 × 256 full logo on white |
| `cfizz-mark-256.png` | 256 × 256 glasses-only mark used at the top of the README |
| `cfizz-social-preview.png` | 1280 × 640 GitHub social preview |

## GitHub social preview

GitHub does not expose a public API for repository social previews. Upload it manually:

1. Open the repository on GitHub.
2. Go to **Settings → Social preview**.
3. Upload `cfizz-social-preview.png`.

The recommended size is 1280 × 640, which this file uses.

## README usage

- `README.md` uses the absolute raw URL for `cfizz-mark-256.png`, so it renders on GitHub and PyPI.
- `README.zh-CN.md` uses the relative path for local preview and GitHub.
