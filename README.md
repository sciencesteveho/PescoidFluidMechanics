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

```python
    # Initialize the class
    dewetvizobj = PescoidVisualizer(
        velocities_file="velocities.mat", area_file="area.csv"
    )

    # Apply smoothing via savgol filter
    dewetvizobj.smooth_area_curve()

    # Apply change point analysis to phase wetting/dewetting
    dewetvizobj.phase_area_curve()

    # Get the value and index where the area of pescoid is at its maximum
    dewetvizobj.max_area
    dewetvizobj.max_area_idx

    # Plot the pescoid's raw area over time
    dewetvizobj.plot_area()

    # Plot the pescoid's smoothed area over time with a vertical line to indicate
    wetting/dewetting phase
    dewetvizobj.plot_area(smoothed=True, phased=True)

    # Plot the rate of the change in pescoid area over time
    dewetvizobj.plot_rate_of_change(smoothed=True)
    dewetvizobj.plot_rate_of_change(smoothed=False)

    # Plot a single line of velocities according to index
    dewetvizobj.plot_single_line_velocity(index=0)

    # Multi-plot velocities for a list of indexes
    dewetvizobj.plot_multiple_line_velocities(indexes=[0, 1, 2])
```


&nbsp;

### Arguments:
| Parameter     |       | Description                           |
|---------------|-------|---------------------------------------|
| area_file | _STR_ | Required |
| velocities_file      | _STR_ | Optional        |
