+++
title = "Milestone 2"
hascode = true
rss = "Description"
rss_title = "Milestone 2"
rss_pubdate = Date(2022, 5, 1)

tags = ["ebm", "solar radiation", "orbital parameters"]
+++

# Milestone 2 - Solar forcing
\toc

So far, we have made the assumption that the solar forcing term remains constant across both space and time. However, this assumption is again too simple, as it overlooks the fact that the Earth's position and orientation in relation to the sun changes throughout the year, leading to a variation in the amount of incoming insolation. 

By incorporating these effects, we can significantly improve the modeling of our solar energy source term. The two primary factors that contribute to the variability of solar radiation received by the Earth are as follows:

1. **Distance from the sun**: Rather than a perfect circle, the Earth's orbit around the sun is an ellipse. Consequently, the distance between the Earth and the sun fluctuates throughout the year, leading to a change in the total amount of solar radiation that reaches the Earth.

2. **Tilted axis of rotation**: Due to the Earth's tilted axis of rotation, different parts of the planet receive varying amounts of solar radiation as it orbits the sun (depending on the latitude). This phenomenon gives rise to the changing seasons.

The derivations presented in this section can be found in the following references:

> [Berger, A. (1978). Long-term variations of daily insolation and Quaternary climatic changes. Journal of Atmospheric Sciences, 35(12), 2362-2367](https://journals.ametsoc.org/downloadpdf/journals/atsc/35/12/1520-0469_1978_035_2362_ltvodi_2_0_co_2.pdf).

> [Berger, A., Loutre, M. F., & Tricot, C. (1993). Insolation and Earth's orbital periods. Journal of Geophysical Research: Atmospheres, 98(D6), 10341-10362](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/93JD00222?casa_token=aTkocQ51OzUAAAAA:dr86J87feYvQiNwbI13BDS8HSHCV4bWn3ouiEkCbl3I9PtuBseVGuK_3YZV8gippC3ZajfOU9wpO-IQ).

### Important definitions

![](/assets/milestone2/Orbits.png)

Before we delve into calculating the impact of the Earth-sun distance and Earth's axis tilting on the solar radiation on Earth, let us define some relevant concepts (see figure above for reference):

* Perihelion: The point in Earth's orbit where it is closest to the Sun. Currently, the perihelion occurs around January 4th.
* Aphelion: The point in Earth's orbit where it is farthest from the Sun. Currently, the aphelion occurs around July 6th.
* Ecliptic: The plane of Earth's orbit around the Sun.
* Northern summer solstice: The point in Earth's orbit where the longest day of the year occurs in the Northern Hemisphere and the shortest day of the year occurs in the Southern Hemisphere. At this point, the plane formed by the rotation axis of Earth and the line that connects Earth and the Sun is perpendicular to the ecliptic. Currently, the Northern summer solstice occurs around June 21st.
* Northern winter solstice: The point in Earth's orbit where the shortest day of the year occurs in the Northern Hemisphere and the longest day of the year occurs in the Southern Hemisphere. At this point, the plane formed by the rotation axis of Earth and the line that connects Earth and the Sun is perpendicular to the ecliptic. Currently, the Northern winter solstice occurs around December 21st.
* Vernal equinox: The point in Earth's orbit where the line that connects Earth and the Sun aligns with the equatorial plane of Earth during the transition from winter to summer in the Northern Hemisphere. This is also the moment when Earth's rotation axis is directly perpendicular to the Sun-Earth line.
* Northern autumnal equinox: The point in Earth's orbit where the line that connects Earth and the Sun aligns with the equatorial plane of Earth during the transition from summer to winter in the Northern Hemisphere. This is also the moment when Earth's rotation axis is directly perpendicular to the Sun-Earth line.
* True longitude of Earth ($\lambda$): Position of Earth at a given time as measured from the vernal equinox. The position $\lambda = 0$ is considered the beginning of an _astronomical year_.

To compute the distance from the sun and Earth's axis tilt angle, we need to consider three parameters that vary over time as a result of the gravitational interactions between Earth and other celestial bodies within the solar system. We call these parameters the **orbital parameters** and define them as

1. Eccentricity ($e$): This parameter determines the degree to which Earth's orbit around the sun deviates from a perfect circle. It is computed as
    $$\label{eq:ecc}
    e = \sqrt{1-\frac{b^2}{a^2}},
    $$ 
    where $a$ represents the semi-major axis, and $b$ represents the semi-minor axis of the ellipse (refer to the figure). Earth's eccentricity undergoes periodic changes over hundreds of thousands of years.
2. Obliquity ($\epsilon$): This parameter, also referred to as axial tilt, represents the angle between Earth's rotational axis and its orbital axis. The obliquity angle is equivalent to the angle between the equatorial plane and the orbital plane. The obliquity angle varies periodically over approximately $\tau_{\epsilon} \approx 20,000$ years.
3.  Precession distance ($\tilde \omega$): Earth's axis rotates with a period of about $40,000$ years, causing the seasons to shift in position within the elliptical orbit of Earth. This phenomenon is known as precession or spin. We define the precession distance as the angle between the aphelion and the vernal equinox.

> **NOTE:** In the paper of Berger, the precession distance ($\tilde \omega$) is sometimes defined from the perihelion and sometimes from the aphelion. We use the aphelion in our implementation as it is the one obtained with equation (6) from Berger.

![](/assets/milestone2/part2-earthspin-nolabel.jpg)
* Figure modified from [https://ugc.berkeley.edu/background-content/earths-spin-tilt-orbit/](https://ugc.berkeley.edu/background-content/earths-spin-tilt-orbit/).

In the following sections, we will assume that the eccentricity, obliquity and precession distance are constant to derive the total solar irradiance.
This is a reasonable assumption for climate simulations in which the orbital parameters are updated at the time interval $\Delta t_{\text{orb}} \ll \tau_{\epsilon}$. 

### Effect of Earth's distance to the sun



<!-- As we saw in the [Radiation Section](#radiation) -->
The inverse-square law states that the total radiation that is received from a source is inversely proportional to the square of the distance from the source.

![](/assets/milestone2/Inverse_square_law.png)
* Source: Wikipedia

As a result, the total radiation that is received at the top layer of Earth's atmosphere is 
$$
S_r(t) = \frac{r_0^2}{r^2} S_0,
$$
where $r$ is the (time-dependent) Earth-sun distance, $r_0 \approx a$ is roughly the mean Earth-sun distance (defined as one anstronomical unit), and $S_0$ is the so-called solar constant, the mean solar elecromagnetic radiation per unit area measured on a surface perpendicular to the solar rays at a distance $r_0$. 

Since Earth's orbit is an ellipse and the sun sits in one focus, we can calulate the Earth-sun distance as
$$\label{eq:r}
r = \frac{a (1-e^2)}{1 - e \cos (\nu)},
$$
where the position (angle) of Earth with respect to the aphelion is
$$\label{eq:nu}
\nu = \lambda - \tilde \omega
$$
The radiation at the top of the atmosphere is then computed as
$$
S_r(t) %= \frac{(1-e \cos (\nu))^2}{(1 - e^2)^2} S_0 
= \frac{(1-e \cos (\lambda - \tilde \omega))^2}{(1 - e^2)^2} S_0.
$$

### Effect of Earth's axial tilt

Due to Earth's axial tilt, the solar irradiance depends on the latitude and the time.

We call the angle between the Earth-sun line and the equatorial plane of Earth _declination angle_ and note it with $\delta$. This angle changes throughout the year, as the Earth revolves around the sun (see figure):

![](/assets/milestone2/Declination_new.png)

In general, the declination angle can be computed as
$$
\sin (\delta) = \sin (\lambda) \sin (\epsilon).
$$

Depending on the declination angle and the latitude, we identify three regions of interest:


#### Latitudes where there is no sunrise (winter)
There is no incoming solar radiation for latitudes that fulfill
$$
|\lat| + |\delta| \ge \frac{\pi}{2} \,\, \text{with} \,\, \lat \delta < 0,
$$
or equivalently
$$
z_0 = -\tan(\lat)\tan(\delta) \ge 1.
$$

For these latitudes, the insolation is simply
$$
S(\lat,t) = 0.
$$


#### Latitudes where there is no sunset (summer)
Some regions of Earth are exposed to incoming solar radiation throughout the entire day during part of the year. These regions fulfill
$$
\label{eq:nosunset}
|\lat| + |\delta| \ge \frac{\pi}{2} \,\, \text{with} \,\, \lat \delta > 0,
$$
or equivalently
$$
z_0 = -\tan(\lat)\tan(\delta) \le -1.
$$

The calculation of the effective insolation is more involved, as it requires the use of spherical astronomy. The result reads as
$$
S(\lat,t) = S_r(t) \sin(\lat) \sin(\delta).
$$

#### Latitudes where there is daily sunrise and sunset
There is daily sunrise and sunset at latitudes that fulfill the condition
$$
- \left( \frac{\pi}{2} - |\delta| \right) < \lat < \frac{\pi}{2} - |\delta|.
$$

For these latitudes, the calculation of the insolation also requires the use of spherical astronomy. The result reads as
$$
S(\lat,t) = S_r(t) \frac{1}{\pi} (H_0 \sin(\lat) \sin(\delta)+\cos(\lat) \cos(\delta) \sin(H_0)),
$$
with the absolute angle of the sun at sunrise and sunset
$$
H_0 = \arccos (z_0).
$$

---

In summary, the insolation is given by the expression
$$
S(\theta,t) =
\begin{cases}
0 & \text{if}~|\lat| + |\delta| \ge \frac{\pi}{2} \,\, \text{with} \,\, \lat \delta < 0,
\\
S_0 \rho(t)\sin(\lat) \sin(\delta) & \text{if}~|\lat| + |\delta| \ge \frac{\pi}{2} \,\, \text{with} \,\, \lat \delta > 0,
 \\
\frac{S_0 \rho(t)}{\pi} (H_0 \sin(\lat) \sin(\delta)+\cos(\lat) \cos(\delta) \sin(H_0)) & \text{if}~- \left( \frac{\pi}{2} - |\delta| \right) < \lat < \frac{\pi}{2} - |\delta|,
\end{cases}
$$
with the normalized distance between Earth and sun squared
$$
\rho(t) =  \left(\frac{1 - e \cos (\nu(t))}{1-e^2}\right)^2.
$$

The total solar forcing is given by
$$
S_{sol}\left(\lat,\varphi,t\right) = S\left(\lat,t\right) (1 - \alpha \left(\lat,\varphi\right)).
$$

### Bonus: Computation of the true longitude $\lambda$

Kepler's second law states that a line segment joining a planet and the Sun sweeps out equal areas during equal intervals of time. Therefore, Earth moves faster when it is closer to the sun. Is a short time interval $\d t$, the area swept is equal to the area of a triangle with base $r \d \nu$ and height $r$. Therefore, we can write Kepler's second law as
$$\label{eq:dAdt}
\frac{\d A}{\d t} = \frac{\pi a b}{T} = \frac{r^2}{2} \frac{\d \nu}{\d t},
$$
where $\pi a b$ is the total area of a ellipse, which is swept in a complete period $T$.

Assuming that the time is measured in years and replacing the definition of the Earth-sun distance \eqref{eq:r} and the eccentricity \eqref{eq:ecc} into \eqref{eq:dAdt} we obtain
$$
\frac{\d \lambda}{\d t} = \frac{2\pi}{(1 - e^2)^{3/2}} (1 - e \cos (\lambda - \tilde \omega))^2.
$$
