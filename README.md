# Pescoid Fluid Dynamics
Tools to plot and analyze wetting and dewetting fluid dynamics in zebrafish early embronic stem cell pescoids.
&nbsp;

## Dependencies
Install via `pip install -r requirements.txt`
```sh
matplotlib==3.5.2
scipy==1.13.0
seaborn==0.11.2
```
&nbsp;

## Usage


&nbsp;

### Required arguments:
| Parameter     |       | Description                           |
|---------------|-------|---------------------------------------|
| --a | _STR_ | Path to fiji trackfile for overlaying |
| --a_name      | _STR_ | Name of fluorescent transgene         |

&nbsp;

### Optional Arguments:
| Parameter           |       | Description                                                 |
|---------------------|-------|-------------------------------------------------------------|
| --peak_detection    |       | Use scipy's peak detection                                  |
| --num_peaks_filter  | _INT_ | Only keep tracks with minimum n number of detected peaks    |
| --periodicity       |       | Plot the periodicity of each trackfile instead of intensity |
| --fourier_transform |       | Plot fourier transformed intensity                          |
| --b       | _STR_ | Path to second fiji trackfile for overlaying                |
| --b_name            | _STR_ | Name of second fluorescent gene for overlaying              |

