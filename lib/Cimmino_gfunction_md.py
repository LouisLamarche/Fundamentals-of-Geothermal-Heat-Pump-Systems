from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.interpolate import interp1d as interp1d
from scipy.constants import pi
from scipy.integrate import quad
from scipy.special import erf,binom
import time as tim

#from mon_pipes import field_thermal_resistance

class Borehole(object):
    """
    Contains information regarding the dimensions and position of a borehole.

    Attributes
    ----------
    H : float
        Borehole length (in meters).
    D : float
        Borehole burried depth (in meters).
    r_b : float
        Borehole radius (in meters).
    x : float
        Position (in meters) of the head of the borehole along the x-axis.
    y : float
        Position (in meters) of the head of the borehole along the y-axis.
    tilt : float
        Angle (in radians) from vertical of the axis of the borehole.
    orientation : float
        Direction (in radians) of the tilt of the borehole.

    """
    def __init__(self, H, D, r_b, x, y, tilt=0., orientation=0.):
        self.H = float(H)      # Borehole length
        self.D = float(D)      # Borehole buried depth
        self.r_b = float(r_b)  # Borehole radius
        self.x = float(x)      # Borehole x coordinate position
        self.y = float(y)      # Borehole y coordinate position
        self.tilt = float(tilt)
        self.orientation = float(orientation)

    def __repr__(self):
        s = ('Borehole(H={self.H}, D={self.D}, r_b={self.r_b}, x={self.x},'
             ' y={self.y}, tilt={self.tilt},'
             ' orientation={self.orientation})').format(self=self)
        return s

    def distance(self, target):
        """
        Evaluate the distance between the current borehole and a target
        borehole.

        Parameters
        ----------
        target : Borehole object
            Target borehole for which the distance is evaluated.

        Returns
        -------
        dis : float
            Distance (in meters) between current borehole and target borehole.

        .. Note::
           The smallest distance returned is equal to the borehole radius.
           This means that the distance between a borehole and itself is
           equal to r_b.

        Examples
        --------
        >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=0., y=0.)
        >>> b2 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
        >>> b1.distance(b2)
        5.0

        """
        dis = max(self.r_b,
                  np.sqrt((self.x - target.x)**2 + (self.y - target.y)**2))
        return dis

    def position(self):
        """
        Returns the position of the borehole.

        Returns
        -------
        pos : tuple
            Position (x, y) (in meters) of the borehole.

        Raises
        ------
        SomeError

        See Also
        --------
        OtherModules

        Examples
        --------
        >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
        >>> b1.position()
        (5.0, 0.0)

        """
        pos = (self.x, self.y)
        return pos


def rectangle_field(N_1, N_2, B_1, B_2, H, D, r_b):
    """
    Build a list of boreholes in a rectangular bore field configuration.

    Parameters
    ----------
    N_1 : int
        Number of borehole in the x direction.
    N_2 : int
        Number of borehole in the y direction.
    B_1 : float
        Distance (in meters) between adjacent boreholes in the x direction.
    B_2 : float
        Distance (in meters) between adjacent boreholes in the y direction.
    H : float
        Borehole length (in meters).
    D : float
        Borehole burried depth (in meters).
    r_b : float
        Borehole radius (in meters).

    Returns
    -------
    boreField : list of Borehole objects
        List of boreholes in the rectangular bore field.

    Examples
    --------
    >>> boreField = gt.boreholes.rectangle_field(N_1=3, N_2=2, B_1=5., B_2=5.,
                                                 H=100., D=2.5, r_b=0.05)

    The bore field is constructed line by line. For N_1=3 and N_2=2, the bore
    field layout is as follows::

     3   4   5

     0   1   2

    """
    borefield = []

    for j in range(N_2):
        for i in range(N_1):
            borefield.append(Borehole(H, D, r_b, x=i*B_1, y=j*B_2))

    return borefield


def L_shaped_field(N_1, N_2, B_1, B_2, H, D, r_b):
    """
    Build a list of boreholes in a L-shaped bore field configuration.

    Parameters
    ----------
    N_1 : int
        Number of borehole in the x direction.
    N_2 : int
        Number of borehole in the y direction.
    B_1 : float
        Distance (in meters) between adjacent boreholes in the x direction.
    B_2 : float
        Distance (in meters) between adjacent boreholes in the y direction.
    H : float
        Borehole length (in meters).
    D : float
        Borehole burried depth (in meters).
    r_b : float
        Borehole radius (in meters).

    Returns
    -------
    boreField : list of Borehole objects
        List of boreholes in the L-shaped bore field.

    Examples
    --------
    >>> boreField = gt.boreholes.L_shaped_field(N_1=3, N_2=2, B_1=5., B_2=5.,
                                                H=100., D=2.5, r_b=0.05)

    The bore field is constructed line by line. For N_1=3 and N_2=2, the bore
    field layout is as follows::

     3

     0   1   2

    """
    borefield = []

    for i in range(N_1):
        borefield.append(Borehole(H, D, r_b, x=i*B_1, y=0.))
    for j in range(1, N_2):
        borefield.append(Borehole(H, D, r_b, x=0., y=j*B_2))

    return borefield


def U_shaped_field(N_1, N_2, B_1, B_2, H, D, r_b):
    """
    Build a list of boreholes in a U-shaped bore field configuration.

    Parameters
    ----------
    N_1 : int
        Number of borehole in the x direction.
    N_2 : int
        Number of borehole in the y direction.
    B_1 : float
        Distance (in meters) between adjacent boreholes in the x direction.
    B_2 : float
        Distance (in meters) between adjacent boreholes in the y direction.
    H : float
        Borehole length (in meters).
    D : float
        Borehole burried depth (in meters).
    r_b : float
        Borehole radius (in meters).

    Returns
    -------
    boreField : list of Borehole objects
        List of boreholes in the U-shaped bore field.

    Examples
    --------
    >>> boreField = gt.boreholes.U_shaped_field(N_1=3, N_2=2, B_1=5., B_2=5.,
                                                H=100., D=2.5, r_b=0.05)

    The bore field is constructed line by line. For N_1=3 and N_2=2, the bore
    field layout is as follows::

     3       4

     0   1   2

    """
    borefield = []

    if N_1 > 2 and N_2 > 1:
        for i in range(N_1):
            borefield.append(Borehole(H, D, r_b, x=i*B_1, y=0.))
        for j in range(1, N_2):
            borefield.append(Borehole(H, D, r_b, x=0, y=j*B_2))
            borefield.append(Borehole(H, D, r_b, x=(N_1-1)*B_1, y=j*B_2))
    else:
        borefield = rectangle_field(N_1, N_2, B_1, B_2, H, D, r_b)

    return borefield


def box_shaped_field(N_1, N_2, B_1, B_2, H, D, r_b):
    """
    Build a list of boreholes in a box-shaped bore field configuration.

    Parameters
    ----------
    N_1 : int
        Number of borehole in the x direction.
    N_2 : int
        Number of borehole in the y direction.
    B_1 : float
        Distance (in meters) between adjacent boreholes in the x direction.
    B_2 : float
        Distance (in meters) between adjacent boreholes in the y direction.
    H : float
        Borehole length (in meters).
    D : float
        Borehole burried depth (in meters).
    r_b : float
        Borehole radius (in meters).

    Returns
    -------
    boreField : list of Borehole objects
        List of boreholes in the box-shaped bore field.

    Examples
    --------
    >>> boreField = gt.boreholes.box_shaped_field(N_1=4, N_2=3, B_1=5., B_2=5.,
                                                  H=100., D=2.5, r_b=0.05)

    The bore field is constructed line by line. For N_1=4 and N_2=3, the bore
    field layout is as follows::

     6   7   8   9

     4           5

     0   1   2   3

    """
    borefield = []


    if N_1 > 2 and N_2 > 2:
        for i in range(N_1):
            borefield.append(Borehole(H, D, r_b, x=i*B_1, y=0.))
        for j in range(1, N_2-1):
            borefield.append(Borehole(H, D, r_b, x=0., y=j*B_2))
            borefield.append(Borehole(H, D, r_b, x=(N_1-1)*B_1, y=j*B_2))
        for i in range(N_1):
            borefield.append(Borehole(H, D, r_b, x=i*B_1, y=(N_2-1)*B_2))
    else:
        borefield = rectangle_field(N_1, N_2, B_1, B_2, H, D, r_b)

    return borefield


def circle_field(N, R, H, D, r_b):
    """
    Build a list of boreholes in a circular field configuration.

    Parameters
    ----------
    N : int
        Number of boreholes in the bore field.
    R : float
        Distance (in meters) of the boreholes from the center of the field.
    H : float
        Borehole length (in meters).
    D : float
        Borehole burried depth (in meters).
    r_b : float
        Borehole radius (in meters).

    Returns
    -------
    boreField : list of Borehole objects
        List of boreholes in the circular shaped bore field.

    Examples
    --------
    >>> boreField = gt.boreholes.circle_field(N=8, R = 5., H=100., D=2.5,
                                              r_b=0.05)

    The bore field is constructed counter-clockwise. For N=8, the bore
    field layout is as follows::

           2
       3       1

     4           0

       5       7
           6

    """
    borefield = []

    for i in range(N):
        borefield.append(Borehole(H, D, r_b, x=R*np.cos(2*pi*i/N),
                                  y=R*np.sin(2*pi*i/N)))

    return borefield


def field_from_file(filename):
    """
    Build a list of boreholes given coordinates and dimensions provided in a
    text file.

    Parameters
    ----------
    filename : str
        Absolute path to text file.

    Returns
    -------
    boreField : list of Borehole objects
        List of boreholes in the bore field.

    The text file should be formatted as follows:

    .. code-block:: none

       # x   y     H     D     r_b
       0.    0.    100.  2.5   0.075
       5.    0.    100.  2.5   0.075
       0.    5.    100.  2.5   0.075
       0.    10.   100.  2.5   0.075
       0.    20.   100.  2.5   0.075

    """
    # Load data from file
    data = np.loadtxt(filename)
    # Build the bore field
    borefield = []
    for line in data:
        x = line[0]
        y = line[1]
        H = line[2]
        D = line[3]
        r_b = line[4]
        borefield.append(Borehole(H, D, r_b, x=x, y=y))

    return borefield


def visualize_field(borefield):
    """
    Plot the top view and 3D view of borehole positions.

    Parameters
    ----------
    borefield : list
        List of boreholes in the bore field.

    Returns
    -------
    fig : figure
        Figure object (matplotlib).

    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    from mpl_toolkits.mplot3d import Axes3D
    # -------------------------------------------------------------------------
    # Initialize figure
    # -------------------------------------------------------------------------
    LW = 1.5    # Line width
    bbox_props = dict(boxstyle="circle,pad=0.3", fc="white", ec="b", lw=LW)

    plt.rc('figure', figsize=(160.0/25.4, 80.0*4.0/4.0/25.4))
    fig = plt.figure()

    # -------------------------------------------------------------------------
    # Top view
    # -------------------------------------------------------------------------
    i = 0   # Initialize borehole index
    ax0 = fig.add_subplot(121)

    for borehole in borefield:
        i += 1  # Increment borehole index
        (x, y) = borehole.position()    # Extract borehole position
        # Add current borehole to the figure
        ax0.plot(x, y, 'k.')
        ax0.text(x, y, i, ha="center", va="center", size=9, bbox=bbox_props)

    # Configure figure axes
    ax0.set_xlabel('x (m)')
    ax0.set_ylabel('y (m)')
    ax0.set_title('Top view')
    plt.axis('equal')
    ax0.xaxis.set_minor_locator(AutoMinorLocator())
    ax0.yaxis.set_minor_locator(AutoMinorLocator())

    # -------------------------------------------------------------------------
    # 3D view
    # -------------------------------------------------------------------------
    i = 0   # Initialize borehole index
    ax1 = fig.add_subplot(122, projection='3d')

    for borehole in borefield:
        i += 1  # Increment borehole index
        # Position of head of borehole
        (x, y) = borehole.position()
        # Position of bottom of borehole
        x_H = x + borehole.H*np.sin(borehole.tilt)*np.cos(borehole.orientation)
        y_H = y + borehole.H*np.sin(borehole.tilt)*np.sin(borehole.orientation)
        z_H = borehole.D + borehole.H*np.cos(borehole.tilt)
        # Add current borehole to the figure
        ax1.plot(np.atleast_1d(x), np.atleast_1d(y), np.atleast_1d(borehole.D),
                 'ko')
        ax1.plot(np.array([x, x_H]),
                 np.array([y, y_H]),
                 -np.array([borehole.D, z_H]), 'k-')

    # Configure figure axes
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_zlabel('z (m)')
    ax1.set_title('3D view')
    plt.axis('equal')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.zaxis.set_minor_locator(AutoMinorLocator())

    plt.tight_layout(rect=[0, 0.0, 0.95, 1.0])

    return fig


def _path_to_inlet(bore_connectivity, bore_index):
    """
    Returns the path from a borehole to the bore field inlet.

    This function raises an error if the supplied borehole connectivity is
    invalid.

    Parameters
    ----------
    bore_connectivity : list
        Index of fluid inlet into each borehole. -1 corresponds to a borehole
        connected to the bore field inlet.
    bore_index : int
        Index of borehole to evaluate path.

    Returns
    -------
    path : list
        List of boreholes leading to the bore field inlet, starting from
        borehole bore_index

    """
    # Initialize path
    path = [bore_index]
    # Index of borehole feeding into borehole (bore_index)
    index_in = bore_connectivity[bore_index]
    # Stop when bore field inlet is reached (index_in == -1)
    while not index_in == -1:
        # Add index of upstream borehole to path
        path.append(index_in)
        # Get index of next upstream borehole
        index_in = bore_connectivity[index_in]

    return path


def _verify_bore_connectivity(bore_connectivity, nBoreholes):
    """
    Verifies that borehole connectivity is valid.

    This function raises an error if the supplied borehole connectivity is
    invalid.

    Parameters
    ----------
    bore_connectivity : list
        Index of fluid inlet into each borehole. -1 corresponds to a borehole
        connected to the bore field inlet.
    nBoreholes : int
        Number of boreholes in the bore field.

    """
    if not len(bore_connectivity) == nBoreholes:
        raise ValueError(
            'The length of the borehole connectivity list does not correspond '
            'to the number of boreholes in the bore field.')
    # Cycle through each borehole and verify that connections lead to -1
    # (-1 is the bore field inlet)
    for i in range(nBoreholes):
        n = 0 # Initialize step counter
        # Index of borehole feeding into borehole i
        index_in = bore_connectivity[i]
        # Stop when bore field inlet is reached (index_in == -1)
        while not index_in == -1:
            index_in = bore_connectivity[index_in]
            n += 1 # Increment step counter
            # Raise error if n exceeds the number of boreholes
            if n > nBoreholes:
                raise ValueError(
                    'The borehole connectivity list is invalid.')
    return


def uniform_heat_extraction(boreholes, time, alpha, use_similarities=True,
                            disTol=0.1, tol=1.0e-6,disp=False):
    """
    Evaluate the g-function with uniform heat extraction along boreholes.

    This function superimposes the finite line source (FLS) solution to
    estimate the g-function of a geothermal bore field.

    Parameters
    ----------
    boreholes : list of Borehole objects
        List of boreholes included in the bore field.
    time : float or array
        Values of time (in seconds) for which the g-function is evaluated.
    alpha : float
        Soil thermal diffusivity (in m2/s).
    use_similarities : bool, optional
        True if similarities are used to limit the number of FLS evaluations.
        Default is True.
    disTol : float, optional
        Absolute tolerance (in meters) on radial distance. Two distances
        (d1, d2) between two pairs of boreholes are considered equal if the
        difference between the two distances (abs(d1-d2)) is below tolerance.
        Default is 0.1.
    tol : float, optional
        Relative tolerance on length and depth. Two lenths H1, H2
        (or depths D1, D2) are considered equal if abs(H1 - H2)/H2 < tol.
        Default is 1.0e-6.

    disp : bool, optional
        Set to true to print progression messages.
        Default is False.

    Returns
    -------
    gFunction : float or array
        Values of the g-function

    Examples
    --------
    >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=0., y=0.)
    >>> b2 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
    >>> alpha = 1.0e-6
    >>> time = np.array([1.0*10**i for i in range(4, 12)])
    >>> gt.gfunction.uniform_heat_extraction([b1, b2], time, alpha)
    array([ 0.75978163,  1.84860837,  2.98861057,  4.33496051,  6.29199383,
        8.13636888,  9.08401497,  9.20736188])

    """
    if disp:
        print(60*'-')
        print('Calculating g-function for uniform heat extraction rate')
        print(60*'-')
    # Initialize chrono
    tic = tim.time()
    # Number of boreholes
    nBoreholes = len(boreholes)
    # Number of time values
    nt = len(np.atleast_1d(time))
    # Initialize heat extraction rates
    Q = np.ones(nBoreholes)
    # Initialize g-function
    gFunction = np.zeros_like(np.atleast_1d(time))
    # Borehole lengths
    H = np.array([b.H for b in boreholes])

    # Calculate borehole to borehole thermal response factors
    h_ij = thermal_response_factors(
        boreholes, time, alpha, use_similarities=use_similarities,
        splitRealAndImage=False, disTol=disTol, tol=tol,disp=disp)
    toc1 = tim.time()

    # Evaluate g-function at all times
    if disp:
        print('Building and solving system of equations ...')
    for i in range(nt):
        Tb = h_ij[:,:,i].dot(Q)
        # The g-function is the average of all borehole wall temperatures
        gFunction[i] = np.dot(Tb, H) / sum(H)
    toc2 = tim.time()

    if disp:
        print('{} sec'.format(toc2 - toc1))
        print('Total time for g-function evaluation: {} sec'.format(
                toc2 - tic))
        print(60*'-')

    # Return float if time is a scalar
    if np.isscalar(time):
        gFunction = np.asscalar(gFunction)

    return gFunction


def uniform_temperature(boreholes, time, alpha, nSegments=12, method='linear',
                        use_similarities=True, disTol=0.1, tol=1.0e-6,disp=False):
    """
    Evaluate the g-function with uniform borehole wall temperature.

    This function superimposes the finite line source (FLS) solution to
    estimate the g-function of a geothermal bore field. Each borehole is
    modeled as a series of finite line source segments, as proposed in
    [#CimminoBernier2014]_.

    Parameters
    ----------
    boreholes : list of Borehole objects
        List of boreholes included in the bore field.
    time : float or array
        Values of time (in seconds) for which the g-function is evaluated.
    alpha : float
        Soil thermal diffusivity (in m2/s).
    nSegments : int, optional
        Number of line segments used per borehole.
        Default is 12.
    method : string, optional
        Interpolation method used for segment-to-segment thermal response
        factors. See documentation for scipy.interpolate.interp1d.
        Default is linear.
    use_similarities : bool, optional
        True if similarities are used to limit the number of FLS evaluations.
        Default is True.
    disTol : float, optional
        Absolute tolerance (in meters) on radial distance. Two distances
        (d1, d2) between two pairs of boreholes are considered equal if the
        difference between the two distances (abs(d1-d2)) is below tolerance.
        Default is 0.1.
    tol : float, optional
        Relative tolerance on length and depth. Two lenths H1, H2
        (or depths D1, D2) are considered equal if abs(H1 - H2)/H2 < tol.
        Default is 1.0e-6.

    disp : bool, optional
        Set to true to print progression messages.
        Default is False.

    Returns
    -------
    gFunction : float or array
        Values of the g-function

    Examples
    --------
    >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=0., y=0.)
    >>> b2 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
    >>> alpha = 1.0e-6
    >>> time = np.array([1.0*10**i for i in range(4, 12)])
    >>> gt.gfunction.uniform_temperature([b1, b2], time, alpha)
    array([ 0.75978079,  1.84859851,  2.98852756,  4.33406497,  6.27830732,
        8.05746656,  8.93697282,  9.04925079])

    References
    ----------
    .. [#CimminoBernier2014] Cimmino, M., & Bernier, M. (2014). A
       semi-analytical method to generate g-functions for geothermal bore
       fields. International Journal of Heat and Mass Transfer, 70, 641-650.

    """
    if disp:
        print(60*'-')
        print('Calculating g-function for uniform borehole wall temperature')
        print(60*'-')
    # Initialize chrono
    tic = tim.time()
    # Number of boreholes
    nBoreholes = len(boreholes)
    # Total number of line sources
    nSources = nSegments*nBoreholes
    # Number of time values
    nt = len(np.atleast_1d(time))
    # Initialize g-function
    gFunction = np.zeros_like(np.atleast_1d(time))
    # Initialize segment heat extraction rates
    Q = np.zeros((nSources, nt))

    # Split boreholes into segments
    boreSegments = _borehole_segments(boreholes, nSegments)
    # Vector of time values
    t = np.atleast_1d(time).flatten()
    # Calculate segment to segment thermal response factors
    h_ij = thermal_response_factors(
        boreSegments, t, alpha, use_similarities=use_similarities,
        splitRealAndImage=True, disTol=disTol, tol=tol,disp=disp)
    toc1 = tim.time()

    if disp:
        print('Building and solving system of equations ...')
    # -------------------------------------------------------------------------
    # Build a system of equation [A]*[X] = [B] for the evaluation of the
    # g-function. [A] is a coefficient matrix, [X] = [Qb,Tb] is a state
    # space vector of the borehole heat extraction rates and borehole wall
    # temperature (equal for all segments), [B] is a coefficient vector.
    # -------------------------------------------------------------------------

    # Segment lengths
    Hb = np.array([b.H for b in boreSegments])
    # Vector of time steps
    dt = np.hstack((t[0], t[1:] - t[:-1]))
    if not np.isscalar(time) and len(time) > 1:
        # Spline object for thermal response factors
        h_dt = interp1d(np.hstack((0., t)),
                        np.dstack((np.zeros((nSources,nSources)), h_ij)),
                        kind=method, axis=2)
        # Thermal response factors evaluated at t=dt
        h_dt = h_dt(dt)
    else:
        h_dt = h_ij
    # Thermal response factor increments
    dh_ij = np.concatenate((h_ij[:,:,0:1], h_ij[:,:,1:]-h_ij[:,:,:-1]), axis=2)

    # Energy conservation: sum([Qb*Hb]) = sum([Hb])
    A_eq2 = np.hstack((Hb, 0.))
    B_eq2 = np.atleast_1d(np.sum(Hb))

    # Build and solve the system of equations at all times
    for p in range(nt):
        # Current thermal response factor matrix
        h_ij_dt = h_dt[:,:,p]
        # Reconstructed load history
        Q_reconstructed = load_history_reconstruction(t[0:p+1], Q[:,0:p+1])
        # Borehole wall temperature for zero heat extraction at current step
        Tb_0 = _temporal_superposition(dh_ij, Q_reconstructed)
        # Spatial superposition: [Tb] = [Tb0] + [h_ij_dt]*[Qb]
        A_eq1 = np.hstack((h_ij_dt, -np.ones((nSources, 1))))
        B_eq1 = -Tb_0
        # Assemble equations
        A = np.vstack((A_eq1, A_eq2))
        B = np.hstack((B_eq1, B_eq2))
        # Solve the system of equations
        X = np.linalg.solve(A, B)
        # Store calculated heat extraction rates
        Q[:,p] = X[0:nSources]
        # The borehole wall temperatures are equal for all segments
        Tb = X[-1]
        gFunction[p] = Tb

    toc2 = tim.time()
    if disp:
        print('{} sec'.format(toc2 - toc1))
        print('Total time for g-function evaluation: {} sec'.format(
                toc2 - tic))
        print(60*'-')

    # Return float if time is a scalar
    if np.isscalar(time):
        gFunction = np.asscalar(gFunction)

    return gFunction,Q


def equal_inlet_temperature(boreholes, UTubes, m_flow, cp, time, alpha,
                            method='linear', nSegments=12,
                            use_similarities=True, disTol=0.1, tol=1.0e-6,
                            disp=False):
    """
    Evaluate the g-function with equal inlet fluid temperatures.

    This function superimposes the finite line source (FLS) solution to
    estimate the g-function of a geothermal bore field. Each borehole is
    modeled as a series of finite line source segments, as proposed in
    [#Cimmino2015]_.

    Parameters
    ----------
    boreholes : list of Borehole objects
        List of boreholes included in the bore field.
    UTubes : list of pipe objects
        Model for pipes inside each borehole.
    m_flow : float or array
        Fluid mass flow rate per borehole (in kg/s).
    cp : float
        Fluid specific isobaric heat capacity (in J/kg.K).
    time : float or array
        Values of time (in seconds) for which the g-function is evaluated.
    alpha : float
        Soil thermal diffusivity (in m2/s).
    nSegments : int, optional
        Number of line segments used per borehole.
        Default is 12.
    method : string, optional
        Interpolation method used for segment-to-segment thermal response
        factors. See documentation for scipy.interpolate.interp1d.
        Default is 'linear'.
    use_similarities : bool, optional
        True if similarities are used to limit the number of FLS evaluations.
        Default is True.
    disTol : float, optional
        Absolute tolerance (in meters) on radial distance. Two distances
        (d1, d2) between two pairs of boreholes are considered equal if the
        difference between the two distances (abs(d1-d2)) is below tolerance.
        Default is 0.1.
    tol : float, optional
        Relative tolerance on length and depth. Two lenths H1, H2
        (or depths D1, D2) are considered equal if abs(H1 - H2)/H2 < tol.
        Default is 1.0e-6.

    disp : bool, optional
        Set to true to print progression messages.
        Default is False.

    Returns
    -------
    gFunction : float or array
        Values of the g-function

    Examples
    --------
    >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=0., y=0.)
    >>> b2 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
    >>> alpha = 1.0e-6
    >>> time = np.array([1.0*10**i for i in range(4, 12)])
    >>> gt.gfunction.uniform_temperature([b1, b2], time, alpha)
    array([ 0.75978079,  1.84859851,  2.98852756,  4.33406497,  6.27830732,
        8.05746656,  8.93697282,  9.04925079])

    References
    ----------
    .. [#Cimmino2015] Cimmino, M. (2015). The effects of borehole thermal
       resistances and fluid flow rate on the g-functions of geothermal bore
       fields. International Journal of Heat and Mass Transfer, 91, 1119-1127.

    """
    if disp:
        print(60*'-')
        print('Calculating g-function for equal inlet fluid temperature')
        print(60*'-')
    # Initialize chrono
    tic = tim.time()
    # Number of boreholes
    nBoreholes = len(boreholes)
    # Total number of line sources
    nSources = nSegments*nBoreholes
    # Number of time values
    nt = len(np.atleast_1d(time))
    # Initialize g-function
    gFunction = np.zeros_like(np.atleast_1d(time))
    # Initialize segment heat extraction rates
    Q = np.zeros((nSources, nt))

    # If m_flow is supplied as float, apply m_flow to all boreholes
    if np.isscalar(m_flow):
        m_flow = np.tile(m_flow, nBoreholes)

    # Split boreholes into segments
    boreSegments = _borehole_segments(boreholes, nSegments)
    # Vector of time values
    t = np.atleast_1d(time).flatten()
    # Calculate segment to segment thermal response factors
    h_ij = thermal_response_factors(
        boreSegments, t, alpha, use_similarities=use_similarities,
        splitRealAndImage=True, disTol=disTol, tol=tol,disp=disp)
    toc1 = tim.time()

    if disp:
        print('Building and solving system of equations ...')
    # -------------------------------------------------------------------------
    # Build a system of equation [A]*[X] = [B] for the evaluation of the
    # g-function. [A] is a coefficient matrix, [X] = [Qb,Tb,Tf_in] is a state
    # space vector of the borehole heat extraction rates, borehole wall
    # temperatures and inlet fluid temperature (equal for all boreholes),
    # [B] is a coefficient vector.
    # -------------------------------------------------------------------------

    # Segment lengths
    Hb = np.array([b.H for b in boreSegments])
    # Vector of time steps
    dt = np.hstack((t[0], t[1:] - t[:-1]))
    if not np.isscalar(time) and len(time) > 1:
        # Spline object for thermal response factors
        h_dt = interp1d(np.hstack((0., t)),
                        np.dstack((np.zeros((nSources,nSources)), h_ij)),
                        kind=method, axis=2)
        # Thermal response factors evaluated at t=dt
        h_dt = h_dt(dt)
    else:
        h_dt = h_ij
    # Thermal response factor increments
    dh_ij = np.concatenate((h_ij[:,:,0:1], h_ij[:,:,1:]-h_ij[:,:,:-1]), axis=2)

    # Energy balance on borehole segments:
    # [Q_{b,i}] = [a_in]*[T_{f,in}] + [a_{b,i}]*[T_{b,i}]
    A_eq2 = np.hstack((-np.eye(nSources), np.zeros((nSources, nSources + 1))))
    B_eq2 = np.zeros(nSources)
    for i in range(nBoreholes):
        # Coefficients for current borehole
        a_in, a_b = UTubes[i].coefficients_borehole_heat_extraction_rate(
                m_flow[i], cp, nSegments)
        # Matrix coefficients ranges
        # Rows
        j1 = i*nSegments
        j2 = (i+1) * nSegments
        # Columns
        n1 = j1 + nSources
        n2 = j2 + nSources
        # Segment length
        Hi = boreholes[i].H / nSegments
        A_eq2[j1:j2, -1:] = a_in / (-2.0*pi*UTubes[i].k_s*Hi)
        A_eq2[j1:j2, n1:n2] = a_b / (-2.0*pi*UTubes[i].k_s*Hi)

    # Energy conservation: sum([Qb*Hb]) = sum([Hb])
    A_eq3 = np.hstack((Hb, np.zeros(nSources + 1)))
    B_eq3 = np.atleast_1d(np.sum(Hb))

    # Build and solve the system of equations at all times
    for p in range(nt):
        # Current thermal response factor matrix
        h_ij_dt = h_dt[:,:,p]
        # Reconstructed load history
        Q_reconstructed = load_history_reconstruction(t[0:p+1], Q[:,0:p+1])
        # Borehole wall temperature for zero heat extraction at current step
        Tb_0 = _temporal_superposition(dh_ij, Q_reconstructed)
        # Spatial superposition: [Tb] = [Tb0] + [h_ij_dt]*[Qb]
        A_eq1 = np.hstack((h_ij_dt,
                           -np.eye(nSources),
                           np.zeros((nSources, 1))))
        B_eq1 = -Tb_0
        # Assemble equations
        A = np.vstack((A_eq1, A_eq2, A_eq3))
        B = np.hstack((B_eq1, B_eq2, B_eq3))
        # Solve the system of equations
        X = np.linalg.solve(A, B)
        # Store calculated heat extraction rates
        Q[:,p] = X[0:nSources]
        # The gFunction is equal to the average borehole wall temperature
        Tb = X[nSources:2*nSources]
        gFunction[p] = Tb.dot(Hb) / np.sum(Hb)

    toc2 = tim.time()
    if disp:
        print('{} sec'.format(toc2 - toc1))
        print('Total time for g-function evaluation: {} sec'.format(
                toc2 - tic))
        print(60*'-')

    # Return float if time is a scalar
    if np.isscalar(time):
        gFunction = np.asscalar(gFunction)

    return gFunction


def mixed_inlet_temperature(boreholes, UTubes, bore_connectivity, m_flow, cp,
                            time, alpha, method='linear', nSegments=12,
                            use_similarities=True, disTol=0.1, tol=1.0e-6,
                            disp=False):
    """
    Evaluate the g-function with mixed inlet fluid temperatures.

    This function superimposes the finite line source (FLS) solution to
    estimate the g-function of a geothermal bore field. Each borehole is
    modeled as a series of finite line source segments, as proposed in
    [#Cimmino2018]_. The piping configurations between boreholes can be any
    combination of series and parallel connections.

    Parameters
    ----------
    boreholes : list of Borehole objects
        List of boreholes included in the bore field.
    UTubes : list of pipe objects
        Model for pipes inside each borehole.
    bore_connectivity : list
        Index of fluid inlet into each borehole. -1 corresponds to a borehole
        connected to the bore field inlet.
    m_flow : float or array
        Fluid mass flow rate in each borehole (in kg/s).
    cp : float
        Fluid specific isobaric heat capacity (in J/kg.K).
    time : float or array
        Values of time (in seconds) for which the g-function is evaluated.
    alpha : float
        Soil thermal diffusivity (in m2/s).
    nSegments : int, optional
        Number of line segments used per borehole.
        Default is 12.
    method : string, optional
        Interpolation method used for segment-to-segment thermal response
        factors. See documentation for scipy.interpolate.interp1d.
        Default is 'linear'.
    use_similarities : bool, optional
        True if similarities are used to limit the number of FLS evaluations.
        Default is True.
    disTol : float, optional
        Absolute tolerance (in meters) on radial distance. Two distances
        (d1, d2) between two pairs of boreholes are considered equal if the
        difference between the two distances (abs(d1-d2)) is below tolerance.
        Default is 0.1.
    tol : float, optional
        Relative tolerance on length and depth. Two lenths H1, H2
        (or depths D1, D2) are considered equal if abs(H1 - H2)/H2 < tol.
        Default is 1.0e-6.

    disp : bool, optional
        Set to true to print progression messages.
        Default is False.

    Returns
    -------
    gFunction : float or array
        Values of the g-function

    Examples
    --------
    >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=0., y=0.)
    >>> b2 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
    >>> Utube1 = gt.pipes.SingleUTube(pos=[(-0.05, 0), (0, -0.05)],
                                      r_in=0.015, r_out=0.02,
                                      borehole=b1,k_s=2, k_g=1, R_fp=0.1)
    >>> Utube2 = gt.pipes.SingleUTube(pos=[(-0.05, 0), (0, -0.05)],
                                      r_in=0.015, r_out=0.02,
                                      borehole=b1,k_s=2, k_g=1, R_fp=0.1)
    >>> bore_connectivity = [-1, 0]
    >>> time = np.array([1.0*10**i for i in range(4, 12)])
    >>> m_flow = 0.25
    >>> cp = 4000.
    >>> alpha = 1.0e-6
    >>> gt.gfunction.mixed_inlet_temperature([b1, b2], [Utube1, Utube2],
                                             bore_connectivity,
                                             m_flow, cp, time, alpha)
    array([ 0.63783569,  1.63305912,  2.72193357,  4.04093857,  5.98242643,
         7.77218495,  8.66198231,  8.77569636])

    References
    ----------
    .. [#Cimmino2018] Cimmino, M. (2018). g-Functions for bore fields with
       mixed parallel and series connections considering the axial fluid
       temperature variations. IGSHPA Research Track, Stockholm. In review.

    """
    if disp:
        print(60*'-')
        print('Calculating g-function for mixed inlet fluid temperatures')
        print(60*'-')
    # Initialize chrono
    tic = tim.time()
    # Number of boreholes
    nBoreholes = len(boreholes)
    # Total number of line sources
    nSources = nSegments*nBoreholes
    # Number of time values
    nt = len(np.atleast_1d(time))
    # Initialize g-function
    gFunction = np.zeros_like(np.atleast_1d(time))
    # Initialize segment heat extraction rates
    Q = np.zeros((nSources, nt))

    # If m_flow is supplied as float, apply m_flow to all boreholes
    if np.isscalar(m_flow):
        m_flow = np.tile(m_flow, nBoreholes)
    m_flow_tot = sum([m_flow[i] for i in range(nBoreholes)
                      if bore_connectivity[i] == -1])

    # Verify that borehole connectivity is valid
    _verify_bore_connectivity(bore_connectivity, nBoreholes)

    # Split boreholes into segments
    boreSegments = _borehole_segments(boreholes, nSegments)
    # Vector of time values
    t = np.atleast_1d(time).flatten()
    # Calculate segment to segment thermal response factors
    h_ij = thermal_response_factors(
        boreSegments, t, alpha, use_similarities=use_similarities,
        splitRealAndImage=True, disTol=disTol, tol=tol,disp=disp)
    toc1 = tim.time()

    if disp:
        print('Building and solving system of equations ...')
    # -------------------------------------------------------------------------
    # g-function. [A] is a coefficient matrix, [X] = [Qb,Tb,Tf_in] is a state
    # Build a system of equation [A]*[X] = [B] for the evaluation of the
    # space vector of the borehole heat extraction rates, borehole wall
    # temperatures and inlet fluid temperature (into the bore field),
    # [B] is a coefficient vector.
    # -------------------------------------------------------------------------

    # Segment lengths
    Hb = np.array([b.H for b in boreSegments])
    # Vector of time steps
    dt = np.hstack((t[0], t[1:] - t[:-1]))
    if not np.isscalar(time) and len(time) > 1:
        # Spline object for thermal response factors
        h_dt = interp1d(np.hstack((0., t)),
                        np.dstack((np.zeros((nSources,nSources)), h_ij)),
                        kind=method, axis=2)
        # Thermal response factors evaluated at t=dt
        h_dt = h_dt(dt)
    else:
        h_dt = h_ij
    # Thermal response factor increments
    dh_ij = np.concatenate((h_ij[:,:,0:1], h_ij[:,:,1:]-h_ij[:,:,:-1]), axis=2)

    # Energy balance on borehole segments:
    # [Q_{b,i}] = [a_in]*[T_{f,in}] + [a_{b,i}]*[T_{b,i}]
    A_eq2 = np.hstack((-np.eye(nSources), np.zeros((nSources, nSources + 1))))
    B_eq2 = np.zeros(nSources)
    for i in range(nBoreholes):
        # Segment length
        Hi = boreholes[i].H / nSegments
        # Rows of equation matrix
        j1 = i*nSegments
        j2 = (i + 1)*nSegments
        # Coefficients for current borehole
        a_in, a_b = UTubes[i].coefficients_borehole_heat_extraction_rate(
                m_flow[i], cp, nSegments)
        # [a_b] is the coefficient matrix for [T_{b,i}]
        n1 = i*nSegments + nSources
        n2 = (i + 1)*nSegments + nSources
        A_eq2[j1:j2, n1:n2] = a_b / (-2.0*pi*UTubes[i].k_s*Hi)

        # Assemble matrix coefficient for [T_{f,in}] and all [T_b]
        path = _path_to_inlet(bore_connectivity, i)
        b_in = a_in
        for j in path[1:]:
            # Coefficients for borehole j
            c_in, c_b = UTubes[j].coefficients_outlet_temperature(
                    m_flow[j], cp, nSegments)
            # Assign the coefficient matrix for [T_{b,j}]
            n2 = (j + 1)*nSegments + nSources
            n1 = j*nSegments + nSources
            A_eq2[j1:j2, n1:n2] = b_in.dot(c_b)/(-2.0*pi*UTubes[i].k_s*Hi)
            # Keep on building coefficient for [T_{f,in}]
            b_in = b_in.dot(c_in)
        A_eq2[j1:j2, -1:] = b_in / (-2.0*pi*UTubes[i].k_s*Hi)

    # Energy conservation: sum([Qb*Hb]) = sum([Hb])
    A_eq3 = np.hstack((Hb, np.zeros(nSources + 1)))
    B_eq3 = np.atleast_1d(np.sum(Hb))

    # Build and solve the system of equations at all times
    for p in range(nt):
        # Current thermal response factor matrix
        h_ij_dt = h_dt[:,:,p]
        # Reconstructed load history
        Q_reconstructed = load_history_reconstruction(t[0:p+1], Q[:,0:p+1])
        # Borehole wall temperature for zero heat extraction at current step
        Tb_0 = _temporal_superposition(dh_ij, Q_reconstructed)
        # Spatial superposition: [Tb] = [Tb0] + [h_ij_dt]*[Qb]
        A_eq1 = np.hstack((h_ij_dt,
                           -np.eye(nSources),
                           np.zeros((nSources, 1))))
        B_eq1 = -Tb_0
        # Assemble equations
        B = np.hstack((B_eq1, B_eq2, B_eq3))
        A = np.vstack((A_eq1, A_eq2, A_eq3))
        # Solve the system of equations
        X = np.linalg.solve(A, B)
        # Store calculated heat extraction rates
        Q[:,p] = X[0:nSources]
        # The gFunction is equal to the average borehole wall temperature
        Tf_in = X[-1]
        Tf_out = Tf_in - 2*pi*UTubes[0].k_s*np.sum(Hb)/(m_flow_tot*cp)
        Tf = 0.5*(Tf_in + Tf_out)
        Rfield = field_thermal_resistance(
                UTubes, bore_connectivity, m_flow, cp)
        Tb_eff = Tf - 2*pi*UTubes[0].k_s*Rfield
        gFunction[p] = Tb_eff

    toc2 = tim.time()

    if disp:
        print('{} sec'.format(toc2 - toc1))
        print('Total time for g-function evaluation: {} sec'.format(
                toc2 - tic))
        print(60*'-')

    # Return float if time is a scalar
    if np.isscalar(time):
        gFunction = np.asscalar(gFunction)

    return gFunction


def load_history_reconstruction(time, Q):
    """
    Reconstructs the load history.

    This function calculates an equivalent load history for an inverted order
    of time step sizes.

    Parameters
    ----------
    time : array
        Values of time (in seconds) in the load history.
    Q : array
        Heat extraction rates (in Watts) of all segments at all times.

    Returns
    -------
    Q_reconstructed : array
        Reconstructed load history.

    """
    # Number of heat sources
    nSources = Q.shape[0]
    # Time step sizes
    dt = np.hstack((time[0], time[1:]-time[:-1]))
    # Time vector
    t = np.hstack((0., time, time[-1] + time[0]))
    # Inverted time step sizes
    dt_reconstructed = dt[::-1]
    # Reconstructed time vector
    t_reconstructed = np.hstack((0., np.cumsum(dt_reconstructed)))
    # Accumulated heat extracted
    f = np.hstack((np.zeros((nSources, 1)), np.cumsum(Q*dt, axis=1)))
    f = np.hstack((f, f[:,-1:]))
    # Create interpolation object for accumulated heat extracted
    sf = interp1d(t, f, kind='linear', axis=1)
    # Reconstructed load history
    Q_reconstructed = (sf(t_reconstructed[1:]) - sf(t_reconstructed[:-1])) \
        / dt_reconstructed

    return Q_reconstructed


def _borehole_segments(boreholes, nSegments):
    """
    Split boreholes into segments.

    This function goes through the list of boreholes and builds a new list,
    with each borehole split into nSegments.

    Parameters
    ----------
    boreholes : list of Borehole objects
        List of boreholes included in the bore field.
    nSegments : int
        Number of line segments used per borehole.

    Returns
    -------
    boreSegments : list
        List of borehole segments.

    """
    boreSegments = []
    for b in boreholes:
        for i in range(nSegments):
            # Divide borehole into segments of equal length
            H = b.H / nSegments
            # Burried depth of the i-th segment
            D = b.D + i * b.H / nSegments
            # Add to list of segments
            boreSegments.append(Borehole(H, D, b.r_b, b.x, b.y))
    return boreSegments


def _temporal_superposition(dh_ij, Q):
    """
    Temporal superposition for inequal time steps.

    Parameters
    ----------
    dh_ij : array
        Values of the segment-to-segment thermal response factor increments at
        the given time step.
    Q_reconstructed : array
        Heat extraction rates of all segments at all times.

    Returns
    -------
    Tb_0 : array
        Current values of borehole wall temperatures assuming no heat
        extraction during current time step.

    """
    # Number of heat sources
    nSources = Q.shape[0]
    # Number of time steps
    nt = Q.shape[1]
    # Borehole wall temperature
    Tb_0 = np.zeros(nSources)
    # Spatial and temporal superpositions
    for it in range(nt):
        Tb_0 += dh_ij[:,:,it].dot(Q[:,nt-it-1])
    return Tb_0


def finite_line_source(
        time, alpha, borehole1, borehole2, reaSource=True, imgSource=True):
    """
    Evaluate the Finite Line Source (FLS) solution.

    This function uses a numerical quadrature to evaluate the one-integral form
    of the FLS solution, as proposed by Claesson and Javed
    [#ClaessonJaved2011]_ and extended to boreholes with different vertical
    positions by Cimmino and Bernier [#CimminoBernier2014]_. The FlS solution
    is given by:

        .. math::
            h_{1\\rightarrow2}(t) &= \\frac{1}{2H_2}
            \\int_{\\frac{1}{\\sqrt{4\\alpha t}}}^{\\infty}
            e^{-d_{12}^2s^2}(I_{real}(s)+I_{imag}(s))ds


            I_{real}(s) &= erfint((D_2-D_1+H_2)s) - erfint((D_2-D_1)s)

            &+ erfint((D_2-D_1-H_1)s) - erfint((D_2-D_1+H_2-H_1)s)

            I_{imag}(s) &= erfint((D_2+D_1+H_2)s) - erfint((D_2+D_1)s)

            &+ erfint((D_2+D_1+H_1)s) - erfint((D_2+D_1+H_2+H_1)s)


            erfint(X) &= \\int_{0}^{X} erf(x) dx

                      &= Xerf(X) - \\frac{1}{\\sqrt{\\pi}}(1-e^{-X^2})

        .. Note::
            The reciprocal thermal response factor
            :math:`h_{2\\rightarrow1}(t)` can be conveniently calculated by:

                .. math::
                    h_{2\\rightarrow1}(t) = \\frac{H_2}{H_1}
                    h_{1\\rightarrow2}(t)

    Parameters
    ----------
    time : float
        Value of time (in seconds) for which the FLS solution is evaluated.
    alpha : float
        Soil thermal diffusivity (in m2/s).
    borehole1 : Borehole object
        Borehole object of the borehole extracting heat.
    borehole2 : Borehole object
        Borehole object for which the FLS is evaluated.
    reaSource : boolean, defaults to True
        True if the real part of the FLS solution is to be included.
    imgSource : boolean, defaults to True
        True if the image part of the FLS solution is to be included.

    Returns
    -------
    h : float
        Value of the FLS solution. The average (over the length) temperature
        drop on the wall of borehole2 due to heat extracted from borehole1 is:

        .. math:: \\Delta T_{b,2} = T_g - \\frac{Q_1}{2\\pi k_s H_2} h

    Examples
    --------
    >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=0., y=0.)
    >>> b2 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
    >>> h = gt.heat_transfer.finite_line_source(4*168*3600., 1.0e-6, b1, b2)
    h = 0.0110473635393

    References
    ----------
    .. [#ClaessonJaved2011] Claesson, J., & Javed, S. (2011). An analytical
       method to calculate borehole fluid temperatures for time-scales from
       minutes to decades. ASHRAE Transactions, 117(2), 279-288.
    .. [#CimminoBernier2014] Cimmino, M., & Bernier, M. (2014). A
       semi-analytical method to generate g-functions for geothermal bore
       fields. International Journal of Heat and Mass Transfer, 70, 641-650.

    """
    def _Ils(s, b1, b2, reaSource, imgSource):
        r = b1.distance(b2)
        func = 0.
        # Function to integrate
        if reaSource:
            # Real part of the FLS solution
            func += _erfint((b2.D - b1.D + b2.H)*s)
            func += -_erfint((b2.D - b1.D)*s)
            func += _erfint((b2.D - b1.D - b1.H)*s)
            func += -_erfint((b2.D - b1.D + b2.H - b1.H)*s)
        if imgSource:
            # Image part of the FLS solution
            func += _erfint((b2.D + b1.D + b2.H)*s)
            func += -_erfint((b2.D + b1.D)*s)
            func += _erfint((b2.D + b1.D + b1.H)*s)
            func += -_erfint((b2.D + b1.D + b2.H + b1.H)*s)
        return 0.5 / (b2.H*s**2) * func * np.exp(-r**2*s**2)

    def _erfint(x):
        # Integral of error function
        return x * erf(x) - 1.0/np.sqrt(np.pi) * (1.0-np.exp(-x**2))

    # Lower bound of integration
    a = 1.0 / np.sqrt(4.0*alpha*time)
    # Evaluate integral using Gauss-Kronrod
    h, err = quad(
        _Ils, a, np.inf, args=(borehole1, borehole2, reaSource, imgSource))
    return h


def thermal_response_factors(
        boreSegments, time, alpha, use_similarities=True,
        splitRealAndImage=True, disTol=0.1, tol=1.0e-6,disp=False):
    """
    Evaluate segment-to-segment thermal response factors.

    This function goes through the list of borehole segments and evaluates
    the segments-to-segment response factors for all times in time.

    Parameters
    ----------
    boreSegments : list of Borehole objects
        List of borehole segments.
    time : float or array
        Values of time (in seconds) for which the g-function is evaluated.
    alpha : float
        Soil thermal diffusivity (in m2/s).
    use_similarities : bool, optional
        True if similarities are used to limit the number of FLS evaluations.
        Default is True.
    splitRealAndImage : bool, optional
        Set to True if similarities are evaluated separately for real and image
        sources. Set to False if similarities are evaluated for the sum of the
        real and image sources.
        Default is True.
    disTol : float, optional
        Absolute tolerance (in meters) on radial distance. Two distances
        (d1, d2) between two pairs of boreholes are considered equal if the
        difference between the two distances (abs(d1-d2)) is below tolerance.
        Default is 0.1.
    tol : float, optional
        Relative tolerance on length and depth. Two lenths H1, H2
        (or depths D1, D2) are considered equal if abs(H1 - H2)/H2 < tol.
        Default is 1.0e-6.

    disp : bool, optional
        Set to true to print progression messages.
        Default is False.

    Returns
    -------
    h_ij : array
        Segment-to-segment thermal response factors.

    """
    # Total number of line sources
    nSources = len(boreSegments)
    # Number of time values
    nt = len(np.atleast_1d(time))

    # Initialize chrono
    tic = tim.time()

    # Initialize segment-to-segment response factors
    h_ij = np.zeros((nSources, nSources, nt))
    # Calculation is based on the choice of use_similarities
    if use_similarities:
        # Calculations with similarities
        if disp:
            print('Identifying similarities ...')
        (nSimPos, simPos, disSimPos, HSimPos, DSimPos,
         nSimNeg, simNeg, disSimNeg, HSimNeg, DSimNeg) = \
            similarities(boreSegments,
                         splitRealAndImage=splitRealAndImage,
                         disTol=disTol,
                         tol=tol)

        toc1 = tim.time()
        if disp:
            print('{} sec'.format(toc1 - tic))
            print('Calculating segment to segment response factors ...')

        # Initialize FLS solution for the real source
        hPos = np.zeros(nt)
        # Initialize FLS solution for the image source
        hNeg = np.zeros(nt)

        # Similarities for real sources
        for s in range(nSimPos):
            n1 = simPos[s][0][0]
            n2 = simPos[s][0][1]
            b1 = boreSegments[n1]
            b2 = boreSegments[n2]
            if splitRealAndImage:
                for it in range(0,nt):
                    ht = finite_line_source(time[it],alpha,b1, b2, reaSource=True, imgSource=False)
                    hPos[it] = ht
            else:
                for it in range(0,nt):
                # FLS solution for combined real and image sources
                    ht = finite_line_source(time[it],alpha,b1, b2, reaSource=True, imgSource=True)
                    hPos[it] = ht
            # Assign thermal response factors to similar segment pairs
            for (i, j) in simPos[s]:
                h_ij[j, i, :] = hPos
                h_ij[i, j, :] = b2.H/b1.H * hPos

        # Similarities for image sources (only if splitRealAndImage=True)
#        if splitRealAndImage:
        for s in range(nSimNeg):
            n1 = simNeg[s][0][0]
            n2 = simNeg[s][0][1]
            b1 = boreSegments[n1]
            b2 = boreSegments[n2]
            # FLS solution for image source only
            for it in range(0,nt):
                ht = finite_line_source(time[it],alpha,b1, b2, reaSource=False, imgSource=True)
                hNeg[it] = ht
            # Assign thermal response factors to similar segment pairs
            for (i, j) in simNeg[s]:
                h_ij[j, i, :] = h_ij[j, i, :] + hNeg
                h_ij[i, j, :] = b2.H/b1.H * h_ij[j, i, :]

    else:
        # Calculations without similarities
        if disp:
            print('Calculating segment to segment response factors ...')
        for i in range(nSources):
            # Segment to same-segment thermal response factor
            # FLS solution for combined real and image sources
            b2 = boreSegments[i]
            h = np.zeros(nt)
            for it in range(0,nt):
                ht = finite_line_source(time[it],alpha,b2,b2)
                h[it] = ht
            # Evaluate the FLS solution at all times in parallel
            h_ij[i, i, :] = h

            # Segment to other segments thermal response factor
            for j in range(i+1, nSources):
                b1 = boreSegments[j]
                # Evaluate the FLS solution at all times in parallel
                h = np.zeros(nt)
                for it in range(0,nt):
                    ht = finite_line_source(time[it],alpha,b1,b2)
                h[it] = ht
                h_ij[i, j, :] = h
                h_ij[j, i, :] = b2.H / b1.H * h_ij[i, j, :]
    toc2 = tim.time()
    if disp:
        print('{} sec'.format(toc2 - tic))

    # Close pool of workers
#    pool.close()
#    pool.join()

    # Return 2d array if time is a scalar
    if np.isscalar(time):
        h_ij = h_ij[:,:,0]

    return h_ij


def similarities(boreholes, splitRealAndImage=True, disTol=0.1, tol=1.0e-6):
    """
    Find similarities in the FLS solution for groups of boreholes.

    This function identifies pairs of boreholes for which the evaluation of the
    Finite Line Source (FLS) solution is equivalent.

    Parameters
    ----------
    boreholes : list of Borehole objects
        List of boreholes for which similarity pairs are identified.
    splitRealAndImage : boolean, defaults to True
        Set to True if similarities are evaluated separately for real and image
        sources. Set to False if similarities are evaluated for the sum of the
        real and image sources.
    distol : float, defaults to 0.1
        Absolute tolerance (in meters) on radial distance. Two distances
        (d1, d2) between two pairs of boreholes are considered equal if the
        difference between the two distances (abs(d1-d2)) is below tolerance.
    tol : float, defaults to 1.0e-6
        Relative tolerance on length and depth. Two lenths H1, H2
        (or depths D1, D2) are considered equal if abs(H1 - H2)/H2 < tol

    Returns
    -------
    nSimPos : integer
        Number of similarities in the evaluation of real sources
        (if splitRealAndImage=True) or sum of real and image sources
        (if splitRealAndImage=False).
    simPos : list of list of tuples
        For each similarity, a list of pairs (tuple) of borehole indexes is
        returned.
    disSimPos : list of floats
        List of distances between boreholes for each similarity.
    HSimPos : list of tuples
        List of lengths of the pairs of boreholes in each similarity.
    DSimPos : list of tuples
        List of depth of the pairs of boreholes in each similarity.
    nSimNeg : integer
        Number of similarities in the evaluation of image sources
        (if splitRealAndImage=True), equals 0 if (splitRealAndImage=False).
    simNeg : list of list of tuples
        For each similarity, a list of pairs (tuple) of borehole indexes is
        returned.
    disSimNeg : list of floats
        List of distances between boreholes for each similarity.
    HSimNeg : list of tuples
        List of lengths of the pairs of boreholes in each similarity.
    DSimNeg : list of tuples
        List of depth of the pairs of boreholes in each similarity.


    Examples
    --------
    >>> b1 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=0., y=0.)
    >>> b2 = gt.boreholes.Borehole(H=150., D=4., r_b=0.075, x=5., y=0.)
    >>> gt.heat_transfer.finite_line_source_similarities([b1, b2])
    2
    [[(0, 0), (1, 1)], [(0, 1)]]
    [0.075, 5.0]
    [(150.0, 150.0), (150.0, 150.0)]
    [(4.0, 4.0), (4.0, 4.0)]
    2
    [[(0, 0), (1, 1)], [(0, 1)]]
    [0.075, 5.0]
    [(150.0, 150.0), (150.0, 150.0)]
    [(4.0, 4.0), (4.0, 4.0)]

    """
    # Initialize pool of workers


    # Group pairs of boreholes by radial distance
    (nDis, disPairs, nPairs, pairs) = \
        _similarities_group_by_distance(boreholes, disTol=disTol)

    # If real and image parts of the FLS are split, evaluate real and image
    # similarities separately:

    if splitRealAndImage:
        realSims = []
        imageSims = []
        # Evaluate similarities for each distance in parallel
        for idis in range(0,nDis):
            rSim  = _similarities_one_distance(pairs[idis],boreholes, kind = 'real')
            realSims.append(rSim)
        # Evaluate similarities for each distance in parallel
            iSim  = _similarities_one_distance(pairs[idis],boreholes,kind = 'image')
            imageSims.append(iSim)

    # Otherwise, evaluate the combined real+image FLS similarities
    else:
        # Evaluate symmetries for each distance in parallel
        realSims = []
        # Evaluate similarities for each distance in parallel
        for idis in range(0,nDis):
            rSim  = _similarities_one_distance(pairs[idis],boreholes, kind = 'real')
            realSims.append(rSim)

    # Close pool of workers

    # Aggregate real similarities for all distances
    nSimPos = 0
    simPos = []
    HSimPos = []
    DSimPos = []
    disSimPos = []
    for i in range(nDis):
        realSim = realSims[i]
        nSim = realSim[0]
        nSimPos += nSim
        disSimPos += [disPairs[i] for j in range(nSim)]
        simPos += realSim[1]
        HSimPos += realSim[2]
        DSimPos += realSim[3]

    # Aggregate image similarities for all distances
    nSimNeg = 0
    simNeg = []
    HSimNeg = []
    DSimNeg = []
    disSimNeg = []
    if splitRealAndImage:
        for i in range(nDis):
            imageSim = imageSims[i]
            nSim = imageSim[0]
            nSimNeg += nSim
            disSimNeg += [disPairs[i] for j in range(nSim)]
            simNeg += imageSim[1]
            HSimNeg += imageSim[2]
            DSimNeg += imageSim[3]

    return nSimPos, simPos, disSimPos, HSimPos, DSimPos, \
            nSimNeg, simNeg, disSimNeg, HSimNeg, DSimNeg


def _similarities_group_by_distance(boreholes, disTol=0.1):
    """
    Groups pairs of boreholes by radial distance between borehole.

    Parameters
    ----------
    boreholes : list of Borehole objects
        List of boreholes in the bore field.
    distol : float, defaults to 0.1
        Absolute tolerance (in meters) on radial distance. Two distances
        (d1, d2) between two pairs of boreholes are considered equal if the
        difference between the two distances (abs(d1-d2)) is below tolerance.

    Returns
    -------
    nDis : int
        Number of unique radial distances between pairs of borehole.
    disPairs : list
        List of radial distances.
    nPairs : list
        List of number of pairs for each radial distance.
    pairs : list
        List of tuples of the borehole indices of borehole pairs at each
        radial distance.

    Raises
    ------
    SomeError

    See Also
    --------
    OtherModules

    Examples
    --------

    """
    # Initialize lists
    nPairs = [1]
    pairs = [[(0, 0)]]
    disPairs = [boreholes[0].r_b]
    nDis = 1

    nb = len(boreholes)
    for i in range(nb):
        b1 = boreholes[i]
        if i == 0:
            i2 = i + 1
        else:
            i2 = i
        for j in range(i2, nb):
            b2 = boreholes[j]
            # Distance between current pair of boreholes
            dis = b1.distance(b2)
            if i == j:
                # The relative tolerance is used for same-borehole
                # distances
                rTol = 1.0e-6 * b1.r_b
            else:
                rTol = disTol
            # Verify if the current pair should be included in the
            # previously identified symmetries
            for k in range(nDis):
                if abs(disPairs[k] - dis) < rTol:
                    pairs[k].append((i, j))
                    nPairs[k] += 1
                    break

            else:
                # Add symmetry to list if no match was found
                nDis += 1
                disPairs.append(dis)
                pairs.append([(i, j)])
                nPairs.append(1)
    return nDis, disPairs, nPairs, pairs


def _similarities_one_distance(pairs, boreholes, kind, tol=1.0e-6):
    """
    Evaluates similarities for all pairs of boreholes separated by the same
    radial distance.

    Parameters
    ----------
    pairs : list
        List of tuples of the borehole indices of borehole pairs at each
        radial distance.
    boreholes : list of Borehole objects
        List of boreholes in the bore field.
    kind : string
        Type of similarity to be evaluated
            - 'real' : similarity in real sources
            - 'image' : similarity in image sources
            - 'realandimage' : similarity for combined real and image sources.
    tol : float, defaults to 1.0e-6
        Relative tolerance on length and depth. Two lenths H1, H2
        (or depths D1, D2) are considered equal if abs(H1 - H2)/H2 < tol

    Returns
    -------
    nSim : int
        Number of similarities.
    sim : list
        For each similarity, a list of pairs (tuple) of borehole indexes is
        returned.
    HSim : list
        List of lengths (tuple) of the pairs of boreholes in each similarity.
    DSim : list
        List of depths (tuple) of the pairs of boreholes in each similarity.

    Raises
    ------
    SomeError

    See Also
    --------
    OtherModules

    Examples
    --------

    """
    # Condition for equivalence of the real part of the FLS solution
    def compare_real_segments(H1a, H1b, H2a, H2b, D1a,
                              D1b, D2a, D2b, tol):
        if (abs((H1a-H1b)/H1a) < tol and
            abs((H2a-H2b)/H2a) < tol and
            abs(((D2a-D1a)-(D2b-D1b))/(D2a-D1a+1e-30)) < tol):
            similarity = True
        else:
            similarity = False
        return similarity

    # Condition for equivalence of the image part of the FLS solution
    def compare_image_segments(H1a, H1b, H2a, H2b,
                               D1a, D1b, D2a, D2b, tol):
        if (abs((H1a-H1b)/H1a) < tol and
            abs((H2a-H2b)/H2a) < tol and
            abs(((D2a+D1a)-(D2b+D1b))/(D2a+D1a+1e-30)) < tol):
            similarity = True
        else:
            similarity = False
        return similarity

    # Condition for equivalence of the full FLS solution
    def compare_realandimage_segments(H1a, H1b, H2a, H2b,
                                      D1a, D1b, D2a, D2b,
                                      tol):
        if (abs((H1a-H1b)/H1a) < tol and
            abs((H2a-H2b)/H2a) < tol and
            abs((D1a-D1b)/(D1a+1e-30)) < tol and
            abs((D2a-D2b)/(D2a+1e-30)) < tol):
            similarity = True
        else:
            similarity = False
        return similarity

    # Initialize comparithermal_response_factorsson function based on input argument
    if kind.lower() == 'real':
        # Check real part of FLS
        compare_segments = compare_real_segments
    elif kind.lower() == 'image':
        # Check image part of FLS
        compare_segments = compare_image_segments
    elif kind.lower() == 'realandimage':
        # Check full real+image FLS
        compare_segments = compare_realandimage_segments

    # Initialize symmetries
    nSim = 1
    pair0 = pairs[0]
    i0 = pair0[0]
    j0 = pair0[1]
    sim = [[pair0]]
    HSim = [(boreholes[i0].H, boreholes[j0].H)]
    DSim = [(boreholes[i0].D, boreholes[j0].D)]

    # Cycle through all pairs of boreholes for the given distance
    for pair in pairs[1:]:
        ibor = pair[0]
        jbor = pair[1]
        b1 = boreholes[ibor]
        b2 = boreholes[jbor]
        # Verify if the current pair should be included in the
        # previously identified symmetries
        for k in range(nSim):
            H1 = HSim[k][0]
            H2 = HSim[k][1]
            D1 = DSim[k][0]
            D2 = DSim[k][1]
            if compare_segments(H1, b1.H, H2, b2.H,
                                D1, b1.D, D2, b2.D, tol):
                sim[k].append((ibor, jbor))
                break
            elif compare_segments(H1, b2.H, H2, b1.H,
                                  D1, b2.D, D2, b1.D, tol):
                sim[k].append((jbor, ibor))
                break

        else:
            # Add symmetry to list if no match was found
            nSim += 1
            sim.append([pair])
            HSim.append((b1.H, b2.H))
            DSim.append((b1.D, b2.D))

    return nSim, sim, HSim, DSim

class _BasePipe(object):
    """
    Template for pipe classes.

    Pipe classes inherit from this class.

    Attributes
    ----------
    borehole : Borehole object
        Borehole class object of the borehole containing the U-Tube.
    nPipes : int
        Number of U-Tubes, equals to 1.
    nInlets : int
        Total number of pipe inlets, equals to 1.
    nOutlets : int
        Total number of pipe outlets, equals to 1.

    """
    def __init__(self, borehole):
        self.b = borehole
        self.nPipes = 1
        self.nInlets = 1
        self.nOutlets = 1

    def get_temperature(self, z, Tin, Tb, m_flow, cp):
        """
        Returns the fluid temperatures of the borehole at a depth (z).

        Parameters
        ----------
        z : float or array
            Depths (in meters) to evaluate the fluid temperatures.
        Tin : float or array
            Inlet fluid temperatures (in Celsius).
        Tb : array
            Borehole wall temperatures (in Celsius).
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).

        Returns
        -------
        Tf : array
            Fluid temperature (in Celsius) in each pipe.

        """
        nSegments = len(np.atleast_1d(Tb))
        z_all = np.atleast_1d(z).flatten()
        Tf = np.zeros((len(z_all), 2*self.nPipes))
        for i in range(len(z_all)):
            zi = z_all[i]
            # Build coefficient matrices
            a_in, a_b = self.coefficients_temperature(zi,
                                                      m_flow,
                                                      cp,
                                                      nSegments)
            # Evaluate fluid temperatures
            Tf[i,:] = a_in.dot(Tin).flatten() + a_b.dot(Tb).flatten()

        # Return 1d array if z was supplied as scalar
        if np.isscalar(z):
            Tf = Tf.flatten()
        return Tf

    def get_inlet_temperature(self, Qf, Tb, m_flow, cp):
        """
        Returns the outlet fluid temperatures of the borehole.

        Parameters
        ----------
        Qf : float or array
            Heat extraction from the fluid circuits (in Watts).
        Tb : float or array
            Borehole wall temperatures (in Celsius).
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).

        Returns
        -------
        Tin : float or array
            Inlet fluid temperatures (in Celsius) into each inlet pipe.

        """
        nSegments = len(np.atleast_1d(Tb))
        # Build coefficient matrices
        a_qf, a_b = self.coefficients_inlet_temperature(m_flow,
                                                        cp,
                                                        nSegments)
        # Evaluate outlet temperatures
        Tin = a_qf.dot(Qf).flatten() + a_b.dot(Tb).flatten()
        # Return float if Tin was supplied as scalar
        if np.isscalar(Qf) and not np.isscalar(Tin):
            Tin = np.asscalar(Tin)
        return Tin

    def get_outlet_temperature(self, Tin, Tb, m_flow, cp):
        """
        Returns the outlet fluid temperatures of the borehole.

        Parameters
        ----------
        Tin : float or array
            Inlet fluid temperatures (in Celsius).
        Tb : float or array
            Borehole wall temperatures (in Celsius).
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).

        Returns
        -------
        Tout : float or array
            Outlet fluid temperatures (in Celsius) from each outlet pipe.

        """
        nSegments = len(np.atleast_1d(Tb))
        # Build coefficient matrices
        a_in, a_b = self.coefficients_outlet_temperature(m_flow,
                                                         cp,
                                                         nSegments)
        # Evaluate outlet temperatures
        Tout = a_in.dot(Tin).flatten() + a_b.dot(Tb).flatten()
        # Return float if Tin was supplied as scalar
        if np.isscalar(Tin) and not np.isscalar(Tout):
            Tout = np.asscalar(Tout)
        return Tout

    def get_borehole_heat_extraction_rate(self, Tin, Tb, m_flow, cp):
        """
        Returns the heat extraction rates of the borehole.

        Parameters
        ----------
        Tin : float or array
            Inlet fluid temperatures (in Celsius).
        Tb : float or array
            Borehole wall temperatures (in Celsius).
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).

        Returns
        -------
        Qb : float or array
            Heat extraction rates along each borehole segment (in Watts).

        """
        nSegments = len(np.atleast_1d(Tb))
        a_in, a_b = self.coefficients_borehole_heat_extraction_rate(m_flow,
                                                                    cp,
                                                                    nSegments)
        Qb = a_in.dot(Tin).flatten() + a_b.dot(Tb).flatten()
        # Return float if Tb was supplied as scalar
        if np.isscalar(Tb) and not np.isscalar(Qb):
            Qb = np.asscalar(Qb)
        return Qb

    def get_fluid_heat_extraction_rate(self, Tin, Tb, m_flow, cp):
        """
        Returns the heat extraction rates of the borehole.

        Parameters
        ----------
        Tin : float or array
            Inlet fluid temperatures (in Celsius).
        Tb : float or array
            Borehole wall temperatures (in Celsius).
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).

        Returns
        -------
        Qf : float or array
            Heat extraction rates from each fluid circuit (in Watts).

        """
        nSegments = len(np.atleast_1d(Tb))
        a_in, a_b = self.coefficients_fluid_heat_extraction_rate(m_flow,
                                                                 cp,
                                                                 nSegments)
        Qf = a_in.dot(Tin).flatten() + a_b.dot(Tb).flatten()
        # Return float if Tb was supplied as scalar
        if np.isscalar(Tin) and not np.isscalar(Qf):
            Qf = np.asscalar(Qf)
        return Qf

    def get_total_heat_extraction_rate(self, Tin, Tb, m_flow, cp):
        """
        Returns the total heat extraction rate of the borehole.

        Parameters
        ----------
        Tin : float or array
            Inlet fluid temperatures (in Celsius).
        Tb : float or array
            Borehole wall temperatures (in Celsius).
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).

        Returns
        -------
        Q : float
            Total net heat extraction rate of the borehole (in Watts).

        """
        Qf = self.get_fluid_heat_extraction_rate(Tin, Tb, m_flow, cp)
        Q = np.sum(Qf)
        return Q

    def coefficients_inlet_temperature(self, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate outlet fluid temperature.

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_{f,in}} = \\mathbf{a_{q,f}} \\mathbf{Q_{f}}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_qf : array
            Array of coefficients for inlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # method_id for coefficients_inlet_temperature is 3
        method_id = 3
        # Check if stored coefficients are available
        if self._check_coefficients(m_flow, cp, nSegments, method_id):
            a_qf, a_b = self._get_stored_coefficients(method_id)
        else:
            # Coefficient matrices for fluid heat extraction rates:
            # [Q_{f}] = [b_in]*[T_{f,in}] + [b_b]*[T_{b}]
            b_in, b_b = self.coefficients_fluid_heat_extraction_rate(m_flow,
                                                                     cp,
                                                                     nSegments)
            b_in_m1 = np.linalg.inv(b_in)

            # Matrices for fluid heat extraction rates:
            # [T_{f,in}] = [a_qf]*[Q_{f}] + [a_b]*[T_{b}]
            a_qf = b_in_m1
            a_b = -b_in_m1.dot(b_b)

            # Store coefficients
            self._set_stored_coefficients(m_flow, cp, nSegments, (a_qf, a_b),
                                          method_id)

        return a_qf, a_b

    def coefficients_outlet_temperature(self, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate outlet fluid temperature.

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_{f,out}} = \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rates (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # method_id for coefficients_outlet_temperature is 4
        method_id = 4
        # Check if stored coefficients are available
        if self._check_coefficients(m_flow, cp, nSegments, method_id):
            a_in, a_b = self._get_stored_coefficients(method_id)
        else:
            # Check if _continuity_condition_base need to be called
            # method_id for _continuity_condition_base is 0
            if self._check_coefficients(m_flow, cp, nSegments, 0):
                b_in, b_out, b_b = self._get_stored_coefficients(0)
            else:
                # Coefficient matrices from continuity condition:
                # [b_out]*[T_{f,out}] = [b_in]*[T_{f,in}] + [b_b]*[T_b]
                b_in, b_out, b_b = self._continuity_condition_base(m_flow,
                                                                   cp,
                                                                   nSegments)

                # Store coefficients
                self._set_stored_coefficients(m_flow, cp, nSegments,
                                              (b_in, b_out, b_b), 0)

            # Final coefficient matrices for outlet temperatures:
            # [T_{f,out}] = [a_in]*[T_{f,in}] + [a_b]*[T_b]
            b_out_m1 = np.linalg.inv(b_out)
            a_in = b_out_m1.dot(b_in)
            a_b = b_out_m1.dot(b_b)

            # Store coefficients
            self._set_stored_coefficients(m_flow, cp, nSegments, (a_in, a_b),
                                          method_id)

        return a_in, a_b

    def coefficients_temperature(self, z, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate fluid temperatures at a depth
        (z).

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_f}(z) = \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        z : float
            Depth (in meters) to evaluate the fluid temperature coefficients.
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # method_id for coefficients_temperature is 5
        method_id = 5

        # Coefficient matrices for outlet temperatures:
        # [T_{f,out}] = [b_in]*[T_{f,in}] + [b_b]*[T_b]
        b_in, b_b = self.coefficients_outlet_temperature(m_flow, cp, nSegments)

        # Check if _continuity_condition_head need to be called
        # method_id for _continuity_condition_head is 1
        if self._check_coefficients(m_flow, cp, nSegments, 1):
            c_in, c_out, c_b = self._get_stored_coefficients(1)
        else:
            # Coefficient matrices for temperatures at depth (z = 0):
            # [T_f](0) = [c_in]*[T_{f,in}] + [c_out]*[T_{f,out}]
            #                              + [c_b]*[T_b]
            c_in, c_out, c_b = self._continuity_condition_head(m_flow,
                                                               cp,
                                                               nSegments)

            # Store coefficients
            self._set_stored_coefficients(m_flow, cp, nSegments,
                                          (c_in, c_out, c_b), 1)

        # Coefficient matrices from general solution:
        # [T_f](z) = [d_f0]*[T_f](0) + [d_b]*[T_b]
        d_f0, d_b = self._general_solution(z, m_flow, cp, nSegments)

        # Final coefficient matrices for temperatures at depth (z):
        # [T_f](z) = [a_in]*[T_{f,in}] + [a_b]*[T_b]
        a_in = d_f0.dot(c_in + c_out.dot(b_in))
        a_b = d_f0.dot(c_b + c_out.dot(b_b)) + d_b

        return a_in, a_b

    def coefficients_borehole_heat_extraction_rate(self,
                                                   m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate heat extraction rates.

        Returns coefficients for the relation:

            .. math::

                \\mathbf{Q_b} = \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # method_id for coefficients_borehole_heat_extraction_rate is 6
        method_id = 6

        nPipes = self.nPipes
        # Check if stored coefficients are available
        if self._check_coefficients(m_flow, cp, nSegments, method_id):
            a_in, a_b = self._get_stored_coefficients(method_id)
        else:
            # Update input variables
            self._format_inputs(m_flow, cp, nSegments)
            m_flow_pipe = self._m_flow_pipe
            cp_pipe = self._cp_pipe
            mcp = np.hstack((-m_flow_pipe[0:nPipes],
                             m_flow_pipe[-nPipes:]))*cp_pipe

            # Initialize coefficient matrices
            a_in = np.zeros((nSegments, self.nInlets))
            a_b = np.zeros((nSegments, nSegments))
            # Heat extraction rates are calculated from an energy balance on a
            # borehole segment.
            z1 = 0.
            aTf1, bTf1 = self.coefficients_temperature(z1,
                                                       m_flow,
                                                       cp,
                                                       nSegments)
            for i in range(nSegments):
                z2 = (i + 1) * self.b.H / nSegments
                aTf2, bTf2 = self.coefficients_temperature(z2,
                                                           m_flow,
                                                           cp,
                                                           nSegments)
                a_in[i, :] = mcp.dot(aTf1 - aTf2)
                a_b[i, :] = mcp.dot(bTf1 - bTf2)
                aTf1, bTf1 = aTf2, bTf2

            # Store coefficients
            self._set_stored_coefficients(m_flow, cp, nSegments, (a_in, a_b),
                                          method_id)

        return a_in, a_b

    def coefficients_fluid_heat_extraction_rate(self, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate heat extraction rates.

        Returns coefficients for the relation:

            .. math::

                \\mathbf{Q_f} = \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # method_id for coefficients_fluid_heat_extraction_rate is 7
        method_id = 7
        # Check if stored coefficients are available
        if self._check_coefficients(m_flow, cp, nSegments, method_id):
            a_in, a_b = self._get_stored_coefficients(method_id)
        else:
            # Update input variables
            self._format_inputs(m_flow, cp, nSegments)

            # Coefficient matrices for outlet temperatures:
            # [T_{f,out}] = [b_in]*[T_{f,in}] + [b_b]*[T_b]
            b_in, b_b = self.coefficients_outlet_temperature(m_flow, cp,
                                                             nSegments)

            # Intermediate matrices for fluid heat extraction rates:
            # [Q_{f}] = [c_in]*[T_{f,in}] + [c_out]*[T_{f,out}]
            MCP = self._m_flow_in * self._cp_in
            c_in = -np.diag(MCP)
            c_out = np.diag(MCP)

            # Matrices for fluid heat extraction rates:
            # [Q_{f}] = [a_in]*[T_{f,in}] + [a_b]*[T_{b}]
            a_in = c_in + c_out.dot(b_in)
            a_b = c_out.dot(b_b)

            # Store coefficients
            self._set_stored_coefficients(m_flow, cp, nSegments, (a_in, a_b),
                                          method_id)

        return a_in, a_b

    def visualize_pipes(self):
        """
        Plot the cross-section view of the borehole.

        Returns
        -------
        fig : figure
            Figure object (matplotlib).

        """
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator

        # Initialize figure
        LW = 1.5    # Line width
        FS = 12.    # Font size

        plt.rc('figure', figsize=(80.0/25.4, 80.0/25.4))
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Borehole wall outline
        borewall = plt.Circle((0., 0.), radius=self.b.r_b,
                              fill=False, linestyle='--', linewidth=LW)
        ax.add_patch(borewall)

        # Pipes
        for i in range(self.nPipes):
            # Coordinates of pipes
            (x_in, y_in) = self.pos[i]
            (x_out, y_out) = self.pos[i + self.nPipes]

            # Pipe outline (inlet)
            pipe_in_in = plt.Circle((x_in, y_in), radius=self.r_in,
                                    fill=False, linestyle='-', linewidth=LW)
            pipe_in_out = plt.Circle((x_in, y_in), radius=self.r_out,
                                     fill=False, linestyle='-', linewidth=LW)
            ax.text(x_in, y_in, i + 1,
                    ha="center", va="center", size=FS)

            # Pipe outline (outlet)
            pipe_out_in = plt.Circle((x_out, y_out), radius=self.r_in,
                                     fill=False, linestyle='-', linewidth=LW)
            pipe_out_out = plt.Circle((x_out, y_out), radius=self.r_out,
                                      fill=False, linestyle='-', linewidth=LW)
            ax.text(x_out, y_out, i + self.nPipes + 1,
                    ha="center", va="center", size=FS)

            ax.add_patch(pipe_in_in)
            ax.add_patch(pipe_in_out)
            ax.add_patch(pipe_out_in)
            ax.add_patch(pipe_out_out)

        # Configure figure axes
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        plt.axis('equal')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        plt.tight_layout()

        return fig

    def _initialize_stored_coefficients(self):
        nMethods = 8    # Number of class methods
        self._stored_coefficients = [() for i in range(nMethods)]
        self._stored_m_flow_cp = [np.empty(self.nInlets)
                                  for i in range(nMethods)]
        self._stored_nSegments = [np.nan for i in range(nMethods)]
        self._m_flow_cp_model_variables = np.empty(self.nInlets)
        self._nSegments_model_variables = np.nan

        return

    def _set_stored_coefficients(self, m_flow, cp, nSegments, coefficients,
                                 method_id):
        self._stored_coefficients[method_id] = coefficients
        self._stored_m_flow_cp[method_id] = m_flow*cp
        self._stored_nSegments[method_id] = nSegments

        return

    def _get_stored_coefficients(self, method_id):
        coefficients = self._stored_coefficients[method_id]

        return coefficients

    def _check_model_variables(self, m_flow, cp, nSegments, tol=1e-6):
        stored_m_flow_cp = self._m_flow_cp_model_variables
        stored_nSegments = self._nSegments_model_variables
        if (np.allclose(m_flow*cp, stored_m_flow_cp, rtol=tol)
                and nSegments == stored_nSegments):
            check = True
        else:
            self._update_model_variables(m_flow, cp, nSegments)
            self._m_flow_cp_model_variables = m_flow*cp
            self._nSegments_model_variables = nSegments
            check = False

        return check

    def _check_coefficients(self, m_flow, cp, nSegments, method_id, tol=1e-6):
        stored_m_flow_cp = self._stored_m_flow_cp[method_id]
        stored_nSegments = self._stored_nSegments[method_id]
        if (np.allclose(m_flow*cp, stored_m_flow_cp, rtol=tol)
                and nSegments == stored_nSegments):
            check = True
        else:
            check = False

        return check

    def _continuity_condition_base(self, m_flow, cp, nSegments):
        """ Returns coefficients for the relation
            [a_out]*[T_{f,out}] = [a_in]*[T_{f,in}] + [a_b]*[T_b]
        """
        raise NotImplementedError(
            '_continuity_condition_base class method not implemented, '
            'this method should return matrices for the relation: '
            '[a_out]*[T_{f,out}] = [a_in]*[T_{f,in}] + [a_b]*[T_b]')

    def _continuity_condition_head(self, m_flow, cp, nSegments):
        """ Returns coefficients for the relation
            [T_f](z=0) = [a_in]*[T_{f,in}] + [a_out]*[T_{f,out}] + [a_b]*[T_b]
        """
        raise NotImplementedError(
            '_continuity_condition_head class method not implemented, '
            'this method should return matrices for the relation: '
            '[T_f](z=0) = [a_in]*[T_{f,in}] + [a_out]*[T_{f,out}] '
            '+ [a_b]*[T_b]')

    def _general_solution(self, z, m_flow, cp, nSegments):
        """ Returns coefficients for the relation
            [T_f](z) = [a_f0]*[T_f](0) + [a_b]*[T_b]
        """
        raise NotImplementedError(
            '_general_solution class method not implemented, '
            'this method should return matrices for the relation: '
            '[T_f](z) = [a_f0]*[T_f](0) + [a_b]*[T_b]')

    def _update_model_variables(self, m_flow, cp, nSegments):
        """
        Evaluate common coefficients needed in other class methods.
        """
        raise NotImplementedError(
            '_update_coefficients class method not implemented, '
            'this method should evaluate common coefficients needed in other '
            'class methods.')

    def _format_inputs(self, m_flow, cp, nSegments):
        """
        Format arrays of mass flow rates and heat capacity.
        """
        raise NotImplementedError(
            '_format_inputs class method not implemented, '
            'this method should format 1d arrays for the inlet mass flow '
            'rates (_m_flow_in), mass flow rates in each pipe (_m_flow_pipe), '
            'heat capacity at each inlet (_cp_in) and heat capacity in each '
            'pipe (_cp_pipe).')


class SingleUTube(_BasePipe):
    """
    Class for single U-Tube boreholes.

    Contains information regarding the physical dimensions and thermal
    characteristics of the pipes and the grout material, as well as methods to
    evaluate fluid temperatures and heat extraction rates based on the work of
    Hellstrom [#Hellstrom1991]_.

    Attributes
    ----------
    pos : list of tuples
        Position (x, y) (in meters) of the pipes inside the borehole.
    r_in : float
        Inner radius (in meters) of the U-Tube pipes.
    r_out : float
        Outter radius (in meters) of the U-Tube pipes.
    borehole : Borehole object
        Borehole class object of the borehole containing the U-Tube.
    k_s : float
        Soil thermal conductivity (in W/m-K).
    k_g : float
        Grout thermal conductivity (in W/m-K).
    R_fp : float
        Fluid to outter pipe wall thermal resistance (m-K/W).
    J : int, optional
        Number of multipoles per pipe to evaluate the thermal resistances.
        Default is 2.
    nPipes : int
        Number of U-Tubes, equals to 1.
    nInlets : int
        Total number of pipe inlets, equals to 1.
    nOutlets : int
        Total number of pipe outlets, equals to 1.

    References
    ----------
    .. [#Hellstrom1991] Hellstrom, G. (1991). Ground heat storage. Thermal
       Analyses of Duct Storage Systems I: Theory. PhD Thesis. University of
       Lund, Department of Mathematical Physics. Lund, Sweden.

    """
    def __init__(self, pos, r_in, r_out, borehole, k_s, k_g, R_fp, J=2):
        self.pos = pos
        self.r_in = r_in
        self.r_out = r_out
        self.b = borehole
        self.k_s = k_s
        self.k_g = k_g
        self.R_fp = R_fp
        self.J = J
        self.nPipes = 1
        self.nInlets = 1
        self.nOutlets = 1

        # Delta-circuit thermal resistances
        self._Rd = thermal_resistances(pos, r_out, borehole.r_b,
                                       k_s, k_g, self.R_fp, J=self.J)[1]

        # Initialize stored_coefficients
        self._initialize_stored_coefficients()

    def _continuity_condition_base(self, m_flow, cp, nSegments):
        """
        Equation that satisfies equal fluid temperatures in both legs of
        each U-tube pipe at depth (z = H).

        Returns coefficients for the relation:

            .. math::

                \\mathbf{a_{out}} T_{f,out} =
                \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_out : array
            Array of coefficients for outlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        # Evaluate coefficient matrices from Hellstrom (1991):
        a_in = ((self._f1(self.b.H) + self._f2(self.b.H))
                / (self._f3(self.b.H) - self._f2(self.b.H)))
        a_in = np.array([[a_in]])

        a_out = np.array([[1.0]])

        a_b = np.zeros((self.nOutlets, nSegments))
        for i in range(nSegments):
            z1 = (nSegments - i - 1) * self.b.H / nSegments
            z2 = (nSegments - i) * self.b.H / nSegments
            dF4 = self._F4(z2) - self._F4(z1)
            dF5 = self._F5(z2) - self._F5(z1)
            a_b[0, i] = (dF4 + dF5) / (self._f3(self.b.H) - self._f2(self.b.H))

        return a_in, a_out, a_b

    def _continuity_condition_head(self, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate fluid temperatures at depth
        (z = 0). These coefficients take into account connections between
        U-tube pipes.

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_f}(z=0) = \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{out}} \\mathbf{T_{f,out}}
                + \\mathbf{a_{b}} \\mathbf{T_{b}}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_out : array
            Array of coefficients for outlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperature.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        # There is only one pipe
        a_in = np.array([[1.0], [0.0]])
        a_out = np.array([[0.0], [1.0]])
        a_b = np.zeros((2, nSegments))

        return a_in, a_out, a_b

    def _general_solution(self, z, m_flow, cp, nSegments):
        """
        General solution for fluid temperatures at a depth (z).

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_f}(z) = \\mathbf{a_{f0}} \\mathbf{T_{f}}(z=0)
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_f0 : array
            Array of coefficients for inlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        a_f0 = np.array([[self._f1(z), self._f2(z)],
                        [-self._f2(z), self._f3(z)]])

        a_b = np.zeros((2*self.nPipes, nSegments))
        N = int(np.ceil(z/self.b.H*nSegments))
        for i in range(N):
            z1 = z - min((i+1)*self.b.H/nSegments, z)
            z2 = z - i * self.b.H / nSegments
            dF4 = self._F4(z2) - self._F4(z1)
            dF5 = self._F5(z2) - self._F5(z1)
            a_b[0, i] = dF4
            a_b[1, i] = -dF5

        return a_f0, a_b

    def _update_model_variables(self, m_flow, cp, nSegments):
        """
        Evaluate dimensionless resistances for Hellstrom (1991) solution.
        """
        # Format mass flow rate and heat capacity inputs
        self._format_inputs(m_flow, cp, nSegments)
        m_flow_in = self._m_flow_in
        cp_in = self._cp_in

        # Dimensionless delta-circuit conductances
        self._beta1 = 1./(self._Rd[0][0]*m_flow_in[0]*cp_in[0])
        self._beta2 = 1./(self._Rd[1][1]*m_flow_in[0]*cp_in[0])
        self._beta12 = 1./(self._Rd[0][1]*m_flow_in[0]*cp_in[0])
        self._beta = 0.5*(self._beta2 - self._beta1)
        # Eigenvalues
        self._gamma = np.sqrt(0.25*(self._beta1+self._beta2)**2
                              + self._beta12*(self._beta1+self._beta2))
        self._delta = 1./self._gamma \
            * (self._beta12 + 0.5*(self._beta1+self._beta2))

    def _format_inputs(self, m_flow, cp, nSegments):
        """
        Format mass flow rate and heat capacity inputs.
        """
        # Format mass flow rate inputs
        if np.isscalar(m_flow):
            # Mass flow rate in each fluid circuit
            m_flow_in = m_flow*np.ones(self.nInlets)
        else:
            # Mass flow rate in each fluid circuit
            m_flow_in = m_flow
        self._m_flow_in = m_flow_in
        # Mass flow rate in pipes
        m_flow_pipe = np.tile(m_flow_in, 2*self.nPipes)
        self._m_flow_pipe = m_flow_pipe

        # Format heat capacity inputs
        if np.isscalar(cp):
            # Heat capacity in each fluid circuit
            cp_in = cp*np.ones(self.nInlets)
        else:
            # Heat capacity in each fluid circuit
            cp_in = cp
        self._cp_in = cp_in
        # Heat capacity in pipes
        cp_pipe = np.tile(cp_in, 2*self.nPipes)
        self._cp_pipe = cp_pipe

    def _f1(self, z):
        """
        Calculate function f1 from Hellstrom (1991)
        """
        f1 = np.exp(self._beta*z)*(np.cosh(self._gamma*z)
                                   - self._delta*np.sinh(self._gamma*z))
        return f1

    def _f2(self, z):
        """
        Calculate function f2 from Hellstrom (1991)
        """
        f2 = np.exp(self._beta*z)*self._beta12/self._gamma \
            * np.sinh(self._gamma*z)
        return f2

    def _f3(self, z):
        """
        Calculate function f3 from Hellstrom (1991)
        """
        f3 = np.exp(self._beta*z)*(np.cosh(self._gamma*z)
                                   + self._delta*np.sinh(self._gamma*z))
        return f3

    def _f4(self, z):
        """
        Calculate function f4 from Hellstrom (1991)
        """
        A = self._delta*self._beta1 + self._beta2*self._beta12/self._gamma
        f4 = np.exp(self._beta*z) \
            * (self._beta1*np.cosh(self._gamma*z) - A*np.sinh(self._gamma*z))
        return f4

    def _f5(self, z):
        """
        Calculate function f5 from Hellstrom (1991)
        """
        B = self._delta*self._beta2 + self._beta1*self._beta12/self._gamma
        f5 = np.exp(self._beta*z) \
            * (self._beta2*np.cosh(self._gamma*z) + B*np.sinh(self._gamma*z))
        return f5

    def _F4(self, z):
        """
        Calculate integral of function f4 from Hellstrom (1991)
        """
        A = self._delta*self._beta1 + self._beta2*self._beta12/self._gamma
        C = self._beta1*self._beta + A*self._gamma
        S = - (self._beta1*self._gamma + self._beta*A)
        denom = (self._beta**2 - self._gamma**2)
        F4 = np.exp(self._beta*z) / denom \
            * (C*np.cosh(self._gamma*z) + S*np.sinh(self._gamma*z))
        return F4

    def _F5(self, z):
        """
        Calculate integral of function f5 from Hellstrom (1991)
        """
        B = self._delta*self._beta2 + self._beta1*self._beta12/self._gamma
        C = self._beta2*self._beta - B*self._gamma
        S = - (self._beta2*self._gamma - self._beta*B)
        denom = (self._beta**2 - self._gamma**2)
        F5 = np.exp(self._beta*z) / denom \
            * (C*np.cosh(self._gamma*z) + S*np.sinh(self._gamma*z))
        return F5


class MultipleUTube(_BasePipe):
    """
    Class for multiple U-Tube boreholes.

    Contains information regarding the physical dimensions and thermal
    characteristics of the pipes and the grout material, as well as methods to
    evaluate fluid temperatures and heat extraction rates based on the work of
    Cimmino [#Cimmino2016]_ for boreholes with any number of U-tubes.

    Attributes
    ----------
    pos : list of tuples
        Position (x, y) (in meters) of the pipes inside the borehole.
    r_in : float
        Inner radius (in meters) of the U-Tube pipes.
    r_out : float
        Outter radius (in meters) of the U-Tube pipes.
    borehole : Borehole object
        Borehole class object of the borehole containing the U-Tube.
    k_s : float
        Soil thermal conductivity (in W/m-K).
    k_g : float
        Grout thermal conductivity (in W/m-K).
    R_fp : float
        Fluid to outter pipe wall thermal resistance (m-K/W).
    J : int, optional
        Number of multipoles per pipe to evaluate the thermal resistances.
        Default is 2.
    nPipes : int
        Number of U-Tubes.
    config : str, defaults to 'parallel'
        Configuration of the U-Tube pipes:
            'parallel' : U-tubes are connected in parallel.
            'series' : U-tubes are connected in series.
    nInlets : int
        Total number of pipe inlets, equals to 1.
    nOutlets : int
        Total number of pipe outlets, equals to 1.

    References
    ----------
    .. [#Cimmino2016] Cimmino, M. (2016). Fluid and borehole wall temperature
       profiles in vertical geothermal boreholes with multiple U-tubes.
       Renewable Energy, 96, 137-147.

    """
    def __init__(self, pos, r_in, r_out, borehole, k_s,
                 k_g, R_fp, nPipes, config='parallel', J=2):
        self.pos = pos
        self.r_in = r_in
        self.r_out = r_out
        self.b = borehole
        self.k_s = k_s
        self.k_g = k_g
        self.R_fp = R_fp
        self.J = J
        self.nPipes = nPipes
        self.nInlets = 1
        self.nOutlets = 1
        self.config = config.lower()

        # Delta-circuit thermal resistances
        self._Rd = thermal_resistances(pos, r_out, borehole.r_b,
                                       k_s, k_g, self.R_fp, J=self.J)[1]

        # Initialize stored_coefficients
        self._initialize_stored_coefficients()

    def _continuity_condition_base(self, m_flow, cp, nSegments):
        """
        Equation that satisfies equal fluid temperatures in both legs of
        each U-tube pipe at depth (z = H).

        Returns coefficients for the relation:

            .. math::

                \\mathbf{a_{out}} T_{f,out} = \\mathbf{a_{in}} T_{f,in}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_out : array
            Array of coefficients for outlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        # Coefficient matrices from continuity condition:
        # [b_u]*[T_{f,u}](z=0) = [b_d]*[T_{f,d}](z=0) + [b_b]*[T_b]
        b_d, b_u, b_b = self._continuity_condition(m_flow, cp, nSegments)
        b_u_m1 = np.linalg.inv(b_u)

        if self.config == 'parallel':
            # Intermediate coefficient matrices:
            # [T_{f,d}](z=0) = [c_in]*[T_{f,in}]
            c_in = np.ones((self.nPipes, 1))

            # Intermediate coefficient matrices:
            # [T_{f,out}] = d_u*[T_{f,u}](z=0)
            mcp = self._m_flow_pipe[-self.nPipes:]*self._cp_pipe[-self.nPipes:]
            d_u = np.reshape(mcp/np.sum(mcp), (1, -1))

            # Final coefficient matrices for continuity at depth (z = H):
            # [a_out][T_{f,out}] = [a_in]*[T_{f,in}] + [a_b]*[T_b]
            a_in = np.linalg.multi_dot([d_u, b_u_m1, b_d, c_in])
            a_out = np.array([[1.0]])
            a_b = np.linalg.multi_dot([d_u, b_u_m1, b_b])

        elif self.config == 'series':
            # Intermediate coefficient matrices:
            # [T_{f,d}](z=0) = [c_in]*[T_{f,in}] + [c_u]*[T_{f,u}](z=0)
            c_in = np.eye(self.nPipes, M=1)
            c_u = np.eye(self.nPipes, k=-1)

            # Intermediate coefficient matrices:
            # [d_u]*[T_{f,u}](z=0) = [d_in]*[T_{f,in}] + [d_b]*[T_b]
            d_u = b_u - b_d.dot(c_u)
            d_in = b_d.dot(c_in)
            d_b = b_b
            d_u_m1 = np.linalg.inv(d_u)

            # Intermediate coefficient matrices:
            # [T_{f,out}] = e_u*[T_{f,u}](z=0)
            e_u = np.eye(self.nPipes, M=1, k=-self.nPipes+1).T

            # Final coefficient matrices for continuity at depth (z = H):
            # [a_out][T_{f,out}] = [a_in]*[T_{f,in}] + [a_b]*[T_b]
            a_in = np.linalg.multi_dot([e_u, d_u_m1, d_in])
            a_out = np.array([[1.0]])
            a_b = np.linalg.multi_dot([e_u, d_u_m1, d_b])

        return a_in, a_out, a_b

    def _continuity_condition_head(self, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate fluid temperatures at depth
        (z = 0). These coefficients take into account connections between
        U-tube pipes.

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_f}(z=0) = \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{out}} \\mathbf{T_{f,out}}
                + \\mathbf{a_{b}} \\mathbf{T_{b}}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_out : array
            Array of coefficients for outlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperature.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        if self.config == 'parallel':
            a_in = np.vstack((np.ones((self.nPipes, self.nInlets)),
                              np.zeros((self.nPipes, self.nInlets))))
            a_out = np.vstack((np.zeros((self.nPipes, self.nOutlets)),
                               np.ones((self.nPipes, self.nOutlets))))
            a_b = np.zeros((2*self.nPipes, nSegments))

        elif self.config == 'series':
            # Coefficient matrices from continuity condition:
            # [b_u]*[T_{f,u}](z=0) = [b_d]*[T_{f,d}](z=0) + [b_b]*[T_b]
            b_d, b_u, b_b = self._continuity_condition(m_flow, cp, nSegments)

            # Intermediate coefficient matrices:
            # [T_{f,d}](z=0) = [c_in]*[T_{f,in}] + [c_u]*[T_{f,u}](z=0)
            c_in = np.eye(self.nPipes, M=1)
            c_u = np.eye(self.nPipes, k=-1)

            # Intermediate coefficient matrices:
            # [d_u]*[T_{f,u}](z=0) = [d_in]*[T_{f,in}] + [d_b]*[T_b]
            d_u = b_u - b_d.dot(c_u)
            d_in = b_d.dot(c_in)
            d_b = b_b
            d_u_m1 = np.linalg.inv(d_u)

            # Intermediate coefficient matrices:
            # [T_f](z=0) = [e_d]*[T_{f,d}](z=0) + [e_u]*[T_{f,u}](z=0)
            e_d = np.eye(2*self.nPipes, M=self.nPipes)
            e_u = np.eye(2*self.nPipes, M=self.nPipes, k=-self.nPipes)

            # Final coefficient matrices for temperatures at depth (z = 0):
            # [T_f](z=0) = [a_in]*[T_{f,in}]+[a_out]*[T_{f,out}]+[a_b]*[T_b]
            a_in = e_d.dot(c_in + np.linalg.multi_dot([c_u, d_u_m1, d_in])) \
                + np.linalg.multi_dot([e_u, d_u_m1, d_in])
            a_out = np.zeros((2*self.nPipes, self.nOutlets))
            a_b = np.linalg.multi_dot([e_d, c_u, d_u_m1, d_b]) \
                + np.linalg.multi_dot([e_u, d_u_m1, d_b])

        return a_in, a_out, a_b

    def _continuity_condition(self, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate fluid temperatures in downward
        and upward flowing pipes at depth (z = 0).

        Returns coefficients for the relation:

            .. math::

                \\mathbf{a_{u}} \\mathbf{T_{f,u}}(z=0) =
                + \\mathbf{a_{d}} \\mathbf{T_{f,d}}(z=0)
                + \\mathbf{a_{b}} \\mathbf{T_{b}}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_d : array
            Array of coefficients for fluid temperature in downward flowing
            pipes.
        a_u : array
            Array of coefficients for fluid temperature in upward flowing
            pipes.
        a_b : array
            Array of coefficients for borehole wall temperature.

        """
        # Load coefficients
        A = self._A
        V = self._V
        Vm1 = self._Vm1
        L = self._L
        Dm1 = self._Dm1

        # Matrix exponential at depth (z = H)
        H = self.b.H
        E = (V.dot(np.diag(np.exp(L*H)))).dot(Vm1)

        # Coefficient matrix for borehole wall temperatures
        IIm1 = np.hstack((np.eye(self.nPipes), -np.eye(self.nPipes)))
        Ones = np.ones((2*self.nPipes, 1))
        a_b = np.zeros((self.nPipes, nSegments))
        for v in range(nSegments):
            z1 = H - v*H/nSegments
            z2 = H - (v + 1)*H/nSegments
            dE = np.diag(np.exp(L*z1) - np.exp(L*z2))
            a_b[:, v:v+1] = np.linalg.multi_dot([IIm1,
                                                 V,
                                                 Dm1,
                                                 dE,
                                                 Vm1,
                                                 A,
                                                 Ones])

        # Configuration-specific inlet and outlet coefficient matrices
        IZER = np.vstack((np.eye(self.nPipes),
                          np.zeros((self.nPipes, self.nPipes))))
        ZERI = np.vstack((np.zeros((self.nPipes, self.nPipes)),
                          np.eye(self.nPipes)))
        a_u = np.linalg.multi_dot([IIm1, E, ZERI])
        a_d = np.linalg.multi_dot([-IIm1, E, IZER])

        return a_d, a_u, a_b

    def _general_solution(self, z, m_flow, cp, nSegments):
        """
        General solution for fluid temperatures at a depth (z).

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_f}(z) = \\mathbf{a_{f0}} \\mathbf{T_{f}}(z=0)
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_f0 : array
            Array of coefficients for inlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        # Load coefficients
        A = self._A
        V = self._V
        Vm1 = self._Vm1
        L = self._L
        Dm1 = self._Dm1

        # Matrix exponential at depth (z)
        a_f0 = (V.dot(np.diag(np.exp(L*z)))).dot(Vm1)

        # Coefficient matrix for borehole wall temperatures
        a_b = np.zeros((2*self.nPipes, nSegments))
        Ones = np.ones((2*self.nPipes, 1))
        for v in range(nSegments):
            dz1 = z - min(z, v*self.b.H/nSegments)
            dz2 = z - min(z, (v + 1)*self.b.H/nSegments)
            E1 = np.diag(np.exp(L*dz1))
            E2 = np.diag(np.exp(L*dz2))
            a_b[:,v:v+1] = np.linalg.multi_dot([V,
                                                Dm1,
                                                E2 - E1,
                                                Vm1,
                                                A,
                                                Ones])

        return a_f0, a_b

    def _update_model_variables(self, m_flow, cp, nSegments):
        """
        Evaluate eigenvalues and eigenvectors for the system of differential
        equations.
        """
        nPipes = self.nPipes
        # Format mass flow rate and heat capacity inputs
        self._format_inputs(m_flow, cp, nSegments)
        m_flow_pipe = self._m_flow_pipe
        cp_pipe = self._cp_pipe

        # Coefficient matrix for differential equations
        self._A = 1.0 / (self._Rd.T * m_flow_pipe * cp_pipe).T
        for i in range(2*nPipes):
            self._A[i, i] = -self._A[i, i] - sum(
                [self._A[i, j] for j in range(2*nPipes) if not i == j])
        for i in range(nPipes, 2*nPipes):
            self._A[i, :] = - self._A[i, :]
        # Eigenvalues and eigenvectors of A
        self._L, self._V = np.linalg.eig(self._A)
        # Inverse of eigenvector matrix
        self._Vm1 = np.linalg.inv(self._V)
        # Diagonal matrix of eigenvalues and inverse
        self._D = np.diag(self._L)
        self._Dm1 = np.diag(1./self._L)

    def _format_inputs(self, m_flow, cp, nSegments):
        """
        Format mass flow rate and heat capacity inputs.
        """
        nPipes = self.nPipes
        # Format mass flow rate inputs
        # Mass flow rate in pipes
        if self.config.lower() == 'parallel':
            m_flow_pipe = np.tile(m_flow/nPipes, 2*self.nPipes)
        elif self.config.lower() == 'series':
            m_flow_pipe = np.tile(m_flow, 2*self.nPipes)
        self._m_flow_pipe = m_flow_pipe
        # Mass flow rate in each fluid circuit
        m_flow_in = np.atleast_1d(m_flow)
        self._m_flow_in = m_flow_in

        # Format heat capacity inputs
        # Heat capacity in each fluid circuit
        cp_in = np.atleast_1d(cp)
        self._cp_in = cp_in
        # Heat capacity in pipes
        cp_pipe = np.tile(cp_in, 2*self.nPipes)
        self._cp_pipe = cp_pipe


class IndependentMultipleUTube(MultipleUTube):
    """
    Class for multiple U-Tube boreholes with independent U-tubes.

    Contains information regarding the physical dimensions and thermal
    characteristics of the pipes and the grout material, as well as methods to
    evaluate fluid temperatures and heat extraction rates based on the work of
    Cimmino [#Cimmino2016b]_ for boreholes with any number of U-tubes.

    Attributes
    ----------
    pos : list of tuples
        Position (x, y) (in meters) of the pipes inside the borehole.
    r_in : float
        Inner radius (in meters) of the U-Tube pipes.
    r_out : float
        Outter radius (in meters) of the U-Tube pipes.
    borehole : Borehole object
        Borehole class object of the borehole containing the U-Tube.
    k_s : float
        Soil thermal conductivity (in W/m-K).
    k_g : float
        Grout thermal conductivity (in W/m-K).
    R_fp : float
        Fluid to outter pipe wall thermal resistance (m-K/W).
    J : int, optional
        Number of multipoles per pipe to evaluate the thermal resistances.
        Default is 2.
    nPipes : int
        Number of U-Tubes.
    nInlets : int
        Total number of pipe inlets, equals to nPipes.
    nOutlets : int
        Total number of pipe outlets, equals to nPipes.

    References
    ----------
    .. [#Cimmino2016b] Cimmino, M. (2016). Fluid and borehole wall temperature
       profiles in vertical geothermal boreholes with multiple U-tubes.
       Renewable Energy, 96, 137-147.

    """
    def __init__(self, pos, r_in, r_out, borehole, k_s,
                 k_g, R_fp, nPipes, J=2):
        self.pos = pos
        self.r_in = r_in
        self.r_out = r_out
        self.b = borehole
        self.k_s = k_s
        self.k_g = k_g
        self.R_fp = R_fp
        self.J = J
        self.nPipes = nPipes
        self.nInlets = nPipes
        self.nOutlets = nPipes

        # Delta-circuit thermal resistances
        self._Rd = thermal_resistances(pos, r_out, borehole.r_b,
                                       k_s, k_g, self.R_fp, J=self.J)[1]

        # Initialize stored_coefficients
        self._initialize_stored_coefficients()

    def _continuity_condition_base(self, m_flow, cp, nSegments):
        """
        Equation that satisfies equal fluid temperatures in both legs of
        each U-tube pipe at depth (z = H).

        Returns coefficients for the relation:

            .. math::

                \\mathbf{a_{out}} T_{f,out} = \\mathbf{a_{in}} T_{f,in}
                + \\mathbf{a_{b}} \\mathbf{T_b}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_out : array
            Array of coefficients for outlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperatures.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        # Coefficient matrices from continuity condition:
        # [b_u]*[T_{f,u}](z=0) = [b_d]*[T_{f,d}](z=0) + [b_b]*[T_b]
        a_in, a_out, a_b = self._continuity_condition(m_flow, cp, nSegments)

        return a_in, a_out, a_b

    def _continuity_condition_head(self, m_flow, cp, nSegments):
        """
        Build coefficient matrices to evaluate fluid temperatures at depth
        (z = 0). These coefficients take into account connections between
        U-tube pipes.

        Returns coefficients for the relation:

            .. math::

                \\mathbf{T_f}(z=0) = \\mathbf{a_{in}} \\mathbf{T_{f,in}}
                + \\mathbf{a_{out}} \\mathbf{T_{f,out}}
                + \\mathbf{a_{b}} \\mathbf{T_{b}}

        Parameters
        ----------
        m_flow : float or array
            Inlet mass flow rate (in kg/s).
        cp : float or array
            Fluid specific isobaric heat capacity (in J/kg.degC).
        nSegments : int
            Number of borehole segments.

        Returns
        -------
        a_in : array
            Array of coefficients for inlet fluid temperature.
        a_out : array
            Array of coefficients for outlet fluid temperature.
        a_b : array
            Array of coefficients for borehole wall temperature.

        """
        # Check if model variables need to be updated
        self._check_model_variables(m_flow, cp, nSegments)

        a_in = np.eye(2*self.nPipes, M=self.nPipes, k=0)
        a_out = np.eye(2*self.nPipes, M=self.nPipes, k=-self.nPipes)
        a_b = np.zeros((2*self.nPipes, nSegments))

        return a_in, a_out, a_b

    def _format_inputs(self, m_flow, cp, nSegments):
        """
        Format mass flow rate and heat capacity inputs.
        """
        # Format mass flow rate inputs
        # Mass flow rate in each fluid circuit
        m_flow_in = np.atleast_1d(m_flow)
        self._m_flow_in = m_flow_in
        # Mass flow rate in pipes
        m_flow_pipe = np.tile(m_flow_in, 2)
        self._m_flow_pipe = m_flow_pipe

        # Format heat capacity inputs
        # Heat capacity in each fluid circuit
        cp_in = np.atleast_1d(cp)
        self._cp_in = cp_in
        # Heat capacity in pipes
        cp_pipe = np.tile(cp_in, 2)
        self._cp_pipe = cp_pipe


def thermal_resistances(pos, r_out, r_b, k_s, k_g, Rfp, J=2):
    """
    Evaluate thermal resistances and delta-circuit thermal resistances.

    This function evaluates the thermal resistances and delta-circuit thermal
    resistances between pipes in a borehole using the multipole method
    [#Claesson2011]_. Thermal resistances are defined by:

    .. math:: \\mathbf{T_f} - T_b = \\mathbf{R} \\cdot \\mathbf{Q_{pipes}}

    Delta-circuit thermal resistances are defined by:

    .. math::

        Q_{i,j} = \\frac{T_{f,i} - T_{f,j}}{R^\\Delta_{i,j}}

        Q_{i,i} = \\frac{T_{f,i} - T_b}{R^\\Delta_{i,i}}

    Parameters
    ----------
    pos : list
        List of positions (x,y) (in meters) of pipes around the center
        of the borehole.
    r_out : float or array
        Outer radius of the pipes (in meters).
    r_b : float
        Borehole radius (in meters).
    k_s : float
        Soil thermal conductivity (in W/m-K).
    k_g : float
        Grout thermal conductivity (in W/m-K).
    Rfp : float or array
        Fluid-to-outer-pipe-wall thermal resistance (in m-K/W).
    J : int, optional
        Number of multipoles per pipe to evaluate the thermal resistances.
        J=1 or J=2 usually gives sufficient accuracy. J=0 corresponds to the
        line source approximation [#Hellstrom1991b]_.
        Default is 2.

    Returns
    -------
    R : array
        Thermal resistances (in m-K/W).
    Rd : array
        Delta-circuit thermal resistances (in m-K/W).

    Examples
    --------
    >>> pos = [(-0.06, 0.), (0.06, 0.)]
    >>> R, Rd = gt.pipes.thermal_resistances(pos, 0.01, 0.075, 2., 1., 0.1,
                                             J=0)
    R = [[ 0.36648149, -0.04855895],
         [-0.04855895,  0.36648149]]
    Rd = [[ 0.31792254, -2.71733044],
          [-2.71733044,  0.31792254]]

    References
    ----------
    .. [#Hellstrom1991b] Hellstrom, G. (1991). Ground heat storage. Thermal
       Analyses of Duct Storage Systems I: Theory. PhD Thesis. University of
       Lund, Department of Mathematical Physics. Lund, Sweden.
    .. [#Claesson2011] Claesson, J., & Hellstrom, G. (2011).
       Multipole method to calculate borehole thermal resistances in a borehole
       heat exchanger. HVAC&R Research, 17(6), 895-911.

    """
    # Number of pipes
    n_p = len(pos)
    # If r_out and/or Rfp are supplied as float, build arrays of size n_p
    if np.isscalar(r_out):
        r_out = np.ones(n_p)*r_out
    if np.isscalar(Rfp):
        Rfp = np.ones(n_p)*Rfp

    R = np.zeros((n_p, n_p))
    if J == 0:
        # Line source approximation
        sigma = (k_g - k_s)/(k_g + k_s)
        for i in range(n_p):
            xi = pos[i][0]
            yi = pos[i][1]
            for j in range(n_p):
                xj = pos[j][0]
                yj = pos[j][1]
                if i == j:
                    # Same-pipe thermal resistance
                    r = np.sqrt(xi**2 + yi**2)
                    R[i, j] = Rfp[i] + 1./(2.*pi*k_g) \
                        *(np.log(r_b/r_out[i]) - sigma*np.log(1 - r**2/r_b**2))
                else:
                    # Pipe to pipe thermal resistance
                    r = np.sqrt((xi-xj)**2 + (yi-yj)**2)
                    ri = np.sqrt(xi**2 + yi**2)
                    rj = np.sqrt(xj**2 + yj**2)
                    dij = np.sqrt((1. - ri**2/r_b**2)*(1.-rj**2/r_b**2) +
                                  r**2/r_b**2)
                    R[i, j] = -1./(2.*pi*k_g) \
                        *(np.log(r/r_b) + sigma*np.log(dij))
    else:
        # Resistances from multipole method are evaluated from the solution of
        # n_p problems
        for m in range(n_p):
            Q_p = np.zeros(n_p)
            Q_p[m] = 1.0
            (T_f, T, it, eps_max) = multipole(pos, r_out, r_b, k_s, k_g,
                                              Rfp, 0., Q_p, J)
            R[:,m] = T_f

    # Delta-circuit thermal resistances
    K = -np.linalg.inv(R)
    for i in range(n_p):
        K[i, i] = -(K[i, i] +
                    sum([K[i, j] for j in range(n_p) if not i == j]))
    Rd = 1.0/K

    return R, Rd


def borehole_thermal_resistance(pipe, m_flow, cp):
    """
    Evaluate the effective borehole thermal resistance.

    Parameters
    ----------
    pipe : pipe object
        Model for pipes inside the borehole.
    m_flow : float
        Fluid mass flow rate (in kg/s).
    cp : float
        Fluid specific isobaric heat capacity (in J/kg.K)

    Returns
    -------
    Rb : float
        Effective borehole thermal resistance (m.K/W).

    """
    # Coefficient for T_{f,out} = a_out*T_{f,in} + [b_out]*[T_b]
    a_out = np.asscalar(
            pipe.coefficients_outlet_temperature(m_flow, cp, nSegments=1)[0])
    # Coefficient for Q_b = [a_Q]*T{f,in} + [b_Q]*[T_b]
    a_Q = np.asscalar(pipe.coefficients_borehole_heat_extraction_rate(
            m_flow, cp, nSegments=1)[0])
    # Borehole length
    H = pipe.b.H
    # Effective borehole thermal resistance
    Rb = -0.5*H*(1. + a_out)/a_Q

    return Rb


def field_thermal_resistance(pipes, bore_connectivity, m_flow, cp):
    """
    Evaluate the effective bore field thermal resistance.

    As proposed in [#Cimmino2018]_.

    Parameters
    ----------
    pipes : list of pipe objects
        Models for pipes inside each borehole.
    bore_connectivity : list
        Index of fluid inlet into each borehole. -1 corresponds to a borehole
        connected to the bore field inlet.
    m_flow : float or array
        Fluid mass flow rate in each borehole (in kg/s).
    cp : float
        Fluid specific isobaric heat capacity (in J/kg.K)

    Returns
    -------
    Rfield : float
        Effective bore field thermal resistance (m.K/W).

    References
    ----------
    .. [#Cimmino2018] Cimmino, M. (2018). g-Functions for bore fields with
       mixed parallel and series connections considering the axial fluid
       temperature variations. IGSHPA Research Track, Stockholm. In review.

    """
    # Number of boreholes
    nBoreholes = len(pipes)
    # If m_flow is supplied as float, apply m_flow to all boreholes
    if np.isscalar(m_flow):
        m_flow = np.tile(m_flow, nBoreholes)
    # Total mass flow rate in the bore field
    m_flow_tot = sum([m_flow[i] for i in range(nBoreholes)
                      if bore_connectivity[i] == -1])

    # Total borehole length
    H_tot = sum([pipes[i].b.H for i in range(nBoreholes)])

    # Verify that borehole connectivity is valid
    _verify_bore_connectivity(bore_connectivity, nBoreholes)


    # Coefficients for T_{f,out} = A_out*T_{f,in} + [B_out]*[T_b], and
    # Q_b = [A_Q]*T{f,in} + [B_Q]*[T_b]
    A_Q = np.array([np.asscalar(
            pipes[i].coefficients_borehole_heat_extraction_rate(
                    m_flow[i], cp, nSegments=1)[0])
            for i in range(len(pipes))])
    for i in range(len(pipes)):
        path = _path_to_inlet(bore_connectivity, i)
        for j in path[1:]:
            a_out = np.asscalar(pipes[j].coefficients_outlet_temperature(
                    m_flow[j], cp, nSegments=1)[0])
            A_Q[i] = A_Q[i]*a_out
    A_out = 1. + np.sum(A_Q)/(m_flow_tot*cp)

    # Effective bore field thermal resistance
    Rfield = -0.5*H_tot*(1. + A_out)/np.sum(A_Q)

    return Rfield


def fluid_friction_factor_circular_pipe(m_flow, r_in, visc, den, epsilon,
                                        tol=1.0e-6):
    """
    Evaluate the Darcy-Weisbach friction factor.

    Parameters
    ----------
    m_flow : float
        Fluid mass flow rate (in kg/s).
    r_in : float
        Inner radius of the pipes (in meters).
    visc : float
        Fluid dynamic viscosity (in kg/m-s).
    den : float
        Fluid density (in kg/m3).
    epsilon : float
        Pipe roughness (in meters).
    tol : float
        Relative convergence tolerance on Darcy friction factor.
        Default is 1.0e-6.

    Returns
    -------
    fDarcy : float
        Darcy friction factor.

    Examples
    --------

    """
    # Hydraulic diameter
    D = 2.*r_in
    # Relative roughness
    E = epsilon / D
    # Fluid velocity
    V_flow = m_flow / den
    A_cs = pi * r_in**2
    V = V_flow / A_cs
    # Reynolds number
    Re = den * V * D / visc

    if Re < 2.3e3:
        # Darcy friction factor for laminar flow
        fDarcy = 64.0 / Re
    else:
        # Colebrook-White equation for rough pipes
        fDarcy = 0.02
        df = 1.0e99
        while abs(df/fDarcy) > tol:
            one_over_sqrt_f = -2.0 * np.log10(E / 3.7
                                              + 2.51/(Re*np.sqrt(fDarcy)))
            fDarcy_new = 1.0 / one_over_sqrt_f**2
            df = fDarcy_new - fDarcy
            fDarcy = fDarcy_new

    return fDarcy


def convective_heat_transfer_coefficient_circular_pipe(m_flow, r_in, visc, den,
                                                       k, cp, epsilon):
    """
    Evaluate the convective heat transfer coefficient for circular pipes.

    Parameters
    ----------
    m_flow : float
        Fluid mass flow rate (in kg/s).
    r_in : float
        Inner radius of the pipes (in meters).
    visc : float
        Fluid dynamic viscosity (in kg/m-s).
    den : float
        Fluid density (in kg/m3).
    k : float
        Fluid thermal conductivity (in W/m-K).
    cp : float
        Fluid specific heat capacity (in J/kg-K).
    epsilon : float
        Pipe roughness (in meters).

    Returns
    -------
    h_fluid : float
        Convective heat transfer coefficient (in W/m2-K).

    Examples
    --------

    """
    # Hydraulic diameter
    D = 2.*r_in
    # Fluid velocity
    V_flow = m_flow / den
    A_cs = pi * r_in**2
    V = V_flow / A_cs
    # Reynolds number
    Re = den * V * D / visc
    # Prandtl number
    Pr = cp * visc / k
    # Darcy friction factor
    fDarcy = fluid_friction_factor_circular_pipe(m_flow, r_in, visc, den,
                                                 epsilon)
    if Re > 2300.:
        # Nusselt number from Gnielinski
        Nu = 0.125*fDarcy * (Re - 1.0e3) * Pr / \
            (1.0 + 12.7 * np.sqrt(0.125*fDarcy) * (Pr**(2.0/3.0) - 1.0))
    else:
        Nu = 3.66
    h_fluid = k * Nu / D

    return h_fluid


def conduction_thermal_resistance_circular_pipe(r_in, r_out, k):
    """
    Evaluate the conduction thermal resistance for circular pipes.

    Parameters
    ----------
    r_in : float
        Inner radius of the pipes (in meters).
    r_out : float
        Outer radius of the pipes (in meters).
    k : float
        Pipe thermal conductivity (in W/m-K).

    Returns
    -------
    R_pipe : float
        Conduction thermal resistance (in m-K/W).

    Examples
    --------

    """
    R_pipe = np.log(r_out/r_in)/(2*pi*k)

    return R_pipe


def multipole(pos, r_p, r_b, k_s, k_g, Rfp, T_b, Q_p, J,
              x_T=np.empty(0), y_T=np.empty(0),
              eps=1e-5, it_max=100):
    """
    Multipole method to calculate borehole thermal resistances in a borehole
    heat exchanger.

    Adapted from the work of Claesson and Hellstrom [#Claesson2011b]_.

    Parameters
    ----------
    pos : list
        List of positions (x,y) (in meters) of pipes around the center
        of the borehole.
    r_p : float or array
        Outer radius of the pipes (in meters).
    r_b : float
        Borehole radius (in meters).
    k_s : float
        Soil thermal conductivity (in W/m-K).
    k_g : float
        Grout thermal conductivity (in W/m-K).
    Rfp : float or array
        Fluid-to-outer-pipe-wall thermal resistance (in m-K/W).
    J : int
        Number of multipoles per pipe to evaluate the thermal resistances.
        J=1 or J=2 usually gives sufficient accuracy. J=0 corresponds to the
        line source approximation.
    Q_p : array
        Thermal energy flows (in W/m) from pipes.
    T_b_av : float
        Average borehole wall temperature (in degC).
    eps : float, optional
        Iteration relative accuracy.
        Default is 1e-5.
    it_max : int, optional
        Maximum number of iterations.
        Default is 100.
    x_T : array, optional
        x-coordinates (in meters) to calculate temperatures.
        Default is np.empty(0).
    y_T : array, optional
        y-coordinates (in meters) to calculate temperatures.
        Default is np.empty(0).

    Returns
    -------
    T_f : array
        Fluid temperatures (in degC) in the pipes.
    T : array
        Requested temperatures (in degC).
    it : int
        Total number of iterations
    eps_max : float
        Maximum error.

    References
    ----------
    .. [#Claesson2011b] Claesson, J., & Hellstrom, G. (2011).
       Multipole method to calculate borehole thermal resistances in a borehole
       heat exchanger. HVAC&R Research, 17(6), 895-911.

    """
    # Pipe coordinates in complex form
    n_p = len(pos)
    z_p = np.array([pos[i][0] + 1.j*pos[i][1] for i in range(n_p)])
    # If r_out and/or Rfp are supplied as float, build arrays of size n_p
    if np.isscalar(r_p):
        r_p = np.ones(n_p)*r_p
    if np.isscalar(Rfp):
        Rfp = np.ones(n_p)*Rfp

    # -------------------------------------
    # Thermal resistance matrix R0 (EQ. 33)
    # -------------------------------------
    pikg = 1.0 / (2.0*pi*k_g)
    sigma = (k_g - k_s)/(k_g + k_s)
    beta_p = 2*pi*k_g*Rfp
    R0 = np.zeros((n_p, n_p))
    for i in range(n_p):
        rbm = r_b**2/(r_b**2 - np.abs(z_p[i])**2)
        R0[i, i] = pikg*(np.log(r_b/r_p[i]) + beta_p[i] + sigma*np.log(rbm))
        for j in range(n_p):
            if i != j:
                dz = np.abs(z_p[i] - z_p[j])
                rbm = r_b**2/np.abs(r_b**2 - z_p[j]*np.conj(z_p[i]))
                R0[i, j] = pikg*(np.log(r_b/dz) + sigma*np.log(rbm))

    # Initialize maximum error and iteration counter
    eps_max = 1.0e99
    it = 0
    # -------------------
    # Multipoles (EQ. 38)
    # -------------------
    if J > 0:
        P = np.zeros((n_p, J), dtype=np.cfloat)
        coeff = -np.array([[(1 - (k+1)*beta_p[m])/(1 + (k+1)*beta_p[m])
                           for k in range(J)] for m in range(n_p)])
        while eps_max > eps and it < it_max:
            it += 1
            eps_max = 0.
            F = _F_mk(Q_p, P, n_p, J, r_b, r_p, z_p, pikg, sigma)
            P_new = coeff*np.conj(F)
            if it == 1:
                diff0 = np.max(np.abs(P_new-P)) - np.min(np.abs(P_new-P))
            diff = np.max(np.abs(P_new-P)) - np.min(np.abs(P_new-P))
            eps_max = diff / diff0
            P = P_new

    # --------------------------
    # Fluid temperatures(EQ. 32)
    # --------------------------
    T_f = T_b + R0.dot(Q_p)
    if J > 0:
        for m in range(n_p):
            dTfm = 0. + 0.j
            for n in range(n_p):
                for j in range(J):
                    # Second term
                    if n != m:
                        dTfm += P[n,j]*(r_p[n]/(z_p[m]-z_p[n]))**(j+1)
                    # Third term
                    dTfm += sigma*P[n,j]*(r_p[n]*np.conj(z_p[m]) \
                                   /(r_b**2 - z_p[n]*np.conj(z_p[m])))**(j+1)
            T_f[m] += np.real(dTfm)

    # -------------------------------
    # Requested temperatures (EQ. 28)
    # -------------------------------
    n_T = len(x_T)
    T = np.zeros(n_T)
    for i in range(n_T):
        z_T = x_T[i] + 1.j*y_T[i]
        dT0 = 0. + 0.j
        dTJ = 0. + 0.j
        for n in range(n_p):
            if np.abs(z_T - z_p[n])/r_p[n] < 1.0:
                # Coordinate inside pipe
                T[i] = T_f[n]
                break
            # Order 0
            if np.abs(z_T) <= r_b:
                # Coordinate inside borehole
                W0 = np.log(r_b/(z_T - z_p[n])) \
                        + sigma*np.log(r_b**2/(r_b**2 - z_p[n]*np.conj(z_T)))
            else:
                # Coordinate outside borehole
                W0 = (1. + sigma)*np.log(r_b/(z_T - z_p[n])) \
                        + sigma*(1. + sigma)/(1. - sigma)*np.log(r_b/z_T)
            dT0 += Q_p[n]*pikg*W0
            # Multipoles
            for j in range(J):
                if np.abs(z_T) <= r_b:
                    # Coordinate inside borehole
                    WJ = (r_p[n]/(z_T - z_p[n]))**(j+1) \
                            + sigma*((r_p[n]*np.conj(z_T))
                                     /(r_b**2 - z_p[n]*np.conj(z_T)))**(j+1)
                else:
                    # Coordinate outside borehole
                    WJ = (1. + sigma)*(r_p[n]/(z_T - z_p[n]))**(j+1)
                dTJ += P[n,j]*WJ
        else:
            T[i] += T_b + np.real(dT0 + dTJ)

    return T_f, T, it, eps_max


def _F_mk(Q_p, P, n_p, J, r_b, r_p, z, pikg, sigma):
    """
    Complex matrix F_mk from Claesson and Hellstrom (2011), EQ. 34
    """
    F = np.zeros((n_p, J), dtype=np.cfloat)
    for m in range(n_p):
        for k in range(J):
            fmk = 0. + 0.j
            for n in range(n_p):
                # First term
                if m != n:
                    fmk += Q_p[n]*pikg/(k+1)*(r_p[m]/(z[n] - z[m]))**(k+1)
                # Second term
                fmk += sigma*Q_p[n]*pikg/(k+1)*(r_p[m]*np.conj(z[n])/(
                        r_b**2 - z[m]*np.conj(z[n])))**(k+1)
                for j in range(J):
                    # Third term
                    if m != n:
                        fmk += P[n,j]*binom(j+k+1, j) \
                                *r_p[n]**(j+1)*(-r_p[m])**(k+1) \
                                /(z[m] - z[n])**(j+k+2)
                    # Fourth term
                    j_pend = np.min((k, j)) + 2
                    for jp in range(j_pend):
                        fmk += sigma*np.conj(P[n,j])*binom(j+1, jp) \
                                *binom(j+k-jp+1, j)*r_p[n]**(j+1) \
                                *r_p[m]**(k+1)*z[m]**(j+1-jp) \
                                *np.conj(z[n])**(k+1-jp) \
                                /(r_b**2 - z[m]*np.conj(z[n]))**(k+j+2-jp)
            F[m,k] = fmk

    return F