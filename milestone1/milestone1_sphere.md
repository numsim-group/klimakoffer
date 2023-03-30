+++
title = "Milestone 1"
hascode = false
rss = "Description"
rss_title = "Milestone 1"
rss_pubdate = Date(2019, 5, 1)

tags = ["climatesystem"]
+++

# Milestone 1 - Spherical coordinates
\toc

The shape of Earth is very close to a sphere with radius of about 6378 km. Therefore, the geometrical model of Earth is reasonably given by the surface of a sphere in our 2D climate model.

There are several conventions for
spherical coordinates. In the geographic coordinate system, we speak of latitude/colatitude and longitude:

* **Latitude** ($\lat$): North-south direction. Latitude lines are parallel to the equator and are assigned the angle from the equator.
* **Colatitude** ($\colat$): Complementary angle from a given latitude (meassured from the north pole).

| Location   | Latitude| Colatitude|
| -----------|---------|-----------|
| North pole | 90°     | 0°        |
| Equator    | 0°      |90°|
| South pole | -90° | 180°|

* **Longitude** ($\varphi$): East-West direction. Longitude lines are perpendicular to the equator with range West -180° to East +180°.

![](/assets/milestone1/LongLat.png)
* **Source:** Wikipedia.

@@colbox-blue
**Remark 1:** The values of longitude and latitude can be given in radians or degrees.
@@

@@colbox-blue
**Remark 2:** In geography, the generation of a map
with long/lat values of a location always
depends on the choice of a reference system: the
so-called **geodetic datum**. The geodetic datum
is a reference ellipsoid. Only in
combination (map+geodetic datum) the
coordinates are precise and can be compared.
@@

We consider the spherical coordinate system, in which all points in the three-dimensional space can be located using three variables:
* Radius ($r \in \R, \, 0 \le r < \infty$),
* Latitude ($\lat \in \R, \, -\pi/2 \le \lat < \pi/2$) **or** colatitude ($\colat \in \R, \, 0 \le \colat < \pi$),
* Longitude: $\varphi \in \R, \, 0 \le \varphi < 2\pi$.

![](/assets/milestone1/SphereCoord.png)

We can map any point in absolute Cartesian coordinates ($x,y,z$) to spherical coordinates ($r,\colat,\varphi$) using the following transformations:
\begin{align}
x &= r \sin \colat \cos \varphi, & r &= \sqrt{x^2 + y^2 + z^2},\\
y &= r \sin \colat \sin \varphi, & \colat &= \arctan \left(\frac{x^2 + y^2}{z} \right),\\
z &= r \cos \colat, & \varphi &= \arctan(y/z).
\end{align}

From mathematical analysis, we know that the Jacobian of the coordinate transformation is given by
\begin{align}
\partialderiv{x}{r} &= \sin \colat \cos \varphi, &
\partialderiv{x}{\varphi} &= -r \sin \colat \sin \varphi, & 
\partialderiv{x}{\colat} &= r \cos \colat \cos \varphi \\
\partialderiv{y}{r} &= \sin \colat \sin \varphi, &
\partialderiv{y}{\varphi} &= r \sin \colat \cos \varphi, & 
\partialderiv{y}{\colat} &= r \cos \colat \sin \varphi \\
\partialderiv{z}{r} &= \sin \colat \cos \colat, &
\partialderiv{z}{\varphi} &= 0, & 
\partialderiv{z}{\colat} &= -r \sin \colat,
\end{align}
respectively in matrix form
\begin{align}
J =
\partialderiv{(x,y,z)}{(r,\colat,\varphi)} =
\begin{bmatrix}
\sin \colat \cos \varphi & -r \sin \colat \sin \varphi & r \cos \colat \cos \varphi \\
\sin \colat \sin \varphi & r \sin \colat \cos \varphi & r \cos \colat \sin \varphi \\
\sin \colat \cos \colat  & 0 & -r \sin \colat 
\end{bmatrix}
\in \R^{3 \times 3},
\end{align}
with the determinant
$$\label{eq:det}
|J|=r^2 \sin \colat.
$$

@@colbox-blue
**Example:**  The coordinate transformation can be used to compute the volume of a sphere with radious R:
\begin{align}
\iiint_{V} \d V &= \iiint_{V} \d x \d y \d z \\
&= \int_{r=0}^R \int_{\colat=0}^{\pi} \int_{\varphi=0}^{2\pi} |J| \d r \d \colat \d \varphi\\
&= \frac{4\pi R^3}{3}
\end{align}
@@
