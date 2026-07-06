# Continuous-gap resolution target for Digital Lensless Holographic Microscopy

This repository contains the MATLAB implementation used for the automated resolution analysis of a continuous-gap resolution target designed for Digital Lensless Holographic Microscopy (DLHM).

The code accompanies the submitted manuscript:

**Continuous-gap resolution target for Digital Lensless Holographic Microscopy**  
Submitted to *Measurement*

**Authors:**  
Maria J. Lopera<sup>1,a,b</sup>, Carlos Buitrago<sup>2,c</sup>, Jorge Sucerquia<sup>3,c</sup>, Yunfeng Nie<sup>4,a</sup>, Heidi Ottevaere<sup>5,a</sup>, and Carlos Trujillo<sup>6,b</sup>

---

## Overview

Resolution assessment in lensless holographic microscopy is commonly performed using discrete targets such as the USAF 1951 resolution chart or radial Siemens-star patterns. While useful, these targets introduce limitations related to discrete spatial-frequency steps, observer-dependent interpretation, and ambiguity in the presence of twin-image artefacts or residual low-frequency background terms.

The continuous-gap target used in this work provides an alternative approach. It consists of pairs of fabricated branches whose edge-to-edge separation varies continuously with position. This allows the resolution limit to be estimated by locating the transition between unresolved and resolved branch profiles, rather than by selecting among discrete resolution groups.

The MATLAB demo included in this repository implements an interactive and semi-automated workflow for estimating the lateral resolution from reconstructed DLHM phase images. The procedure includes background removal, user-tunable contrast enhancement, Canny edge detection, Hough-transform-based center estimation, circular cross-section extraction, adjacent peak-pair analysis, and conversion of the detected resolution radius into a physical resolution value.

---

## Method summary

Starting from a reconstructed phase image of the target, the algorithm first removes low-frequency background contributions and enhances the visibility of the fabricated branches. Edge-based processing is then used to detect the dominant target geometry and estimate the target center. If the automatic center detection is not reliable, the user can manually select the center from the graphical interface.

Once the center is defined, the image is sampled along a sequence of concentric circular trajectories. Each circular cross-section produces a one-dimensional angular profile. Intersections between the circular trajectory and the fabricated branches appear as local peaks in this profile. For each radius, the algorithm identifies adjacent peak pairs and evaluates whether the two branches can be considered resolved according to peak prominence, separation, valley depth, and local profile shape.

The first radius at which two adjacent peaks are reliably distinguishable is interpreted as the resolution threshold. This radius is converted into a target coordinate using the effective object-plane pixel size,

\[
p_x = \frac{p_\mathrm{sensor}}{M},
\]

where \(p_\mathrm{sensor}\) is the sensor pixel size and \(M\) is the geometrical magnification of the DLHM system. The corresponding coordinate is then mapped to the calibrated continuous-gap equation,

\[
s(x)=Ax^3+Bx^2+Cx+10^{-6},
\]

with

\[
A = 290682.4 \ \mathrm{m^{-3}}, \quad
B = -100.5 \ \mathrm{m^{-2}}, \quad
C = 0.0305 \ \mathrm{m^{-1}}.
\]

The resulting value \(s(x)\) is reported as the measured lateral resolution in micrometers.

---

## Repository purpose

This repository is intended to provide a transparent and reproducible implementation of the resolution-estimation procedure described in the manuscript. The code is organized as a demo-style analysis tool, with intermediate plots included to help users inspect each step of the processing pipeline.

The graphical interface allows the user to:

- tune the background-removal and contrast-enhancement parameters,
- inspect the Canny edge map and detected Hough lines,
- estimate the target center automatically or select it manually,
- define the maximum circular cross-section radius,
- visualize individual angular cross-sections,
- scan all radii and identify resolved adjacent peak pairs,
- compute the final resolution in micrometers,
- export plots and analysis results.

Example input images, metadata files, and target information are provided to reproduce the analysis shown in the manuscript.
