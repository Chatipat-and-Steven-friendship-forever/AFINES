"""Read and write AFINES output files."""

import numpy as np


# -----------------------------------------------------------------------------
# [load/save]
# load and save AFINES output files


def load_pe(dirname):
    """
    Load box parameters, potential energies, and virials.

    Parameters
    ----------
    dirname : string
        AFINES output directory.

    Returns
    -------
    data : (N,) ndarray of pe_dtype
        Box parameters, potential energies, and virials at each frame.

    """
    return np.loadtxt("{}/data/pe.txt".format(dirname), dtype=pe_dtype)


def save_pe(dirname, data):
    """
    Save box parameters, potential energies, and virials.

    Parameters
    ----------
    dirname : string
        AFINES output directory.
    data : (N,) ndarray of pe_dtype
        Box parameters, potential energies, and virials at each frame.

    """
    np.savetxt("{}/data/pe.txt".format(dirname), data, delimiter="\t")


def load_filament_e(dirname):
    """
    Load per filament potential energies.

    Parameters
    ----------
    dirname : string
        AFINES output directory.

    Returns
    -------
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of filament_e_dtype
        Per filament potential energies.

    """
    return load_file(
        "{}/data/filament_e.txt".format(dirname), filament_e_dtype
    )


def save_filament_e(dirname, times, frames):
    """
    Save per filament potential energies.

    Parameters
    ----------
    dirname : string
        AFINES output directory.
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of filament_e_dtype
        Per filament potential energies.

    """
    save_file("{}/data/filament_e.txt".format(dirname), times, frames)


def load_actins(dirname):
    """
    Load actin coordinates.

    Parameters
    ----------
    dirname : string
        AFINES output directory.

    Returns
    -------
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of actins_dtype
        Actin coordinates.

    """
    return load_file("{}/txt_stack/actins.txt".format(dirname), actins_dtype)


def save_actins(dirname, times, frames):
    """
    Save actin coordinates.

    Parameters
    ----------
    dirname : string
        AFINES output directory.
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of actins_dtype
        Actin coordinates.

    """
    save_file("{}/txt_stack/actins.txt".format(dirname), times, frames)


def load_links(dirname):
    """
    Load link coordinates.

    Parameters
    ----------
    dirname : string
        AFINES output directory.

    Returns
    -------
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of links_dtype
        Link coordinates.

    """
    return load_file("{}/txt_stack/links.txt".format(dirname), links_dtype)


def save_links(dirname, times, frames):
    """
    Save link coordinates.

    Parameters
    ----------
    dirname : string
        AFINES output directory.
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of links_dtype
        Link coordinates.

    """
    save_file("{}/txt_stack/links.txt".format(dirname), times, frames)


def load_amotors(dirname):
    """
    Load motor coordinates and attachment locations.

    Parameters
    ----------
    dirname : string
        AFINES output directory.

    Returns
    -------
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of motors_dtype
        Motor coordinates and attachment locations.

    """
    return load_file("{}/txt_stack/amotors.txt".format(dirname), motors_dtype)


def save_amotors(dirname, times, frames):
    """
    Save motor coordinates and attachment locations.

    Parameters
    ----------
    dirname : string
        AFINES output directory.
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of motors_dtype
        Motor coordinates and attachment locations.

    """
    save_file("{}/txt_stack/amotors.txt".format(dirname), times, frames)


def load_pmotors(dirname):
    """
    Load crosslinker coordinates and attachment locations.

    Parameters
    ----------
    dirname : string
        AFINES output directory.

    Returns
    -------
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of motors_dtype
        Crosslinker coordinates and attachment locations.

    """
    return load_file("{}/txt_stack/pmotors.txt".format(dirname), motors_dtype)


def save_pmotors(dirname, times, frames):
    """
    Save crosslinker coordinates and attachment locations.

    Parameters
    ----------
    dirname : string
        AFINES output directory.
    times : (N,) ndarray of float
        Time of each frame.
    frames : (N, M) ndarray of motors_dtype
        Crosslinker coordinates and attachment locations.

    """
    save_file("{}/txt_stack/pmotors.txt".format(dirname), times, frames)


def load_file(filename, dtype):
    """
    Load an AFINES output file.

    An AFINES output file consists of a sequence of frames.
    Each frame starts with a single line specifying
    the time of the frame and the number of records in the frame,
    formatted as

    ``t = ``*float*``\tN =``*int*

    This is followed by the records,
    each formatted as a single line of tab separated values.

    Parameters
    ----------
    filename : string
        Path to AFINES output file.
    dtype : dtype
        Data type of records.

    Returns
    -------
    times : float
        Time of each frame.
    frames : (N, M) ndarray of dtype
        List of frames. Each frame is a list of records.

    """
    times = []
    frames = []
    with open(filename) as f:
        for line in f:
            # frame header
            line = line.split()
            if not (
                len(line) == 6
                and line[0] == "t"
                and line[1] == "="
                and line[3] == "N"
                and line[4] == "="
            ):
                raise ValueError("frame header format is invalid")
            time = float(line[2])
            n = int(line[5])

            # frame records
            frame = np.loadtxt(f, dtype=dtype, max_rows=n)
            if len(frame) != n:
                raise ValueError("frame at time {} is incomplete".format(time))

            times.append(time)
            frames.append(frame)
    return np.array(times), np.array(frames)


def save_file(filename, times, frames):
    """
    Save an AFINES output file.

    An AFINES output file consists of a sequence of frames.
    Each frame starts with a single line specifying
    the time of the frame and the number of records in the frame,
    formatted as

    ``t = ``*float*``\tN =``*int*

    This is followed by the records,
    each formatted as a single line of tab separated values.

    Parameters
    ----------
    filename : string
        Path to AFINES output file.
    times : float
        Time of each frame.
    frames : (N, M) ndarray of dtype
        List of frames. Each frame is a list of records.

    """
    if len(times) != len(frames):
        raise ValueError(
            "number of times ({}) and frames ({}) don't match".format(
                len(times), len(frames)
            )
        )
    with open(filename, "w") as f:
        for time, frame in zip(times, frames):
            f.write("t = {}\tN = {}\n".format(time, len(frame)))
            np.savetxt(f, frame, delimiter="\t")


# -----------------------------------------------------------------------------
# [convert]
# convert between numpy arrays with the output dtype
# and lists of filament/motor coordinates


def actins2coords(actins):
    """
    Convert actin frames to coordinates.

    Parameters
    ----------
    actins : (N M,) ndarray of actins_dtype
        List of frames. Each frame is a list of records.

    Returns
    -------
    coords : (N, M, 2) ndarray of float
        List of filaments. Each filament is a list of actin coordinates.

    """
    coords = []
    last_index = None
    for actin in actins:
        if actin["f_index"] != last_index:
            coords.append([])
            last_index = actin["f_index"]
        coords[-1].append((actin["x"], actin["y"]))
    coords = [np.array(filament) for filament in coords]
    return np.array(coords)


def coords2actins(coords, radius):
    """
    Convert actin coordinates to frames.

    Parameters
    ----------
    coords : (N, M, 2) ndarray of float
        List of filaments. Each filament is a list of actin coordinates.

    radius : float
        Actin radius.

    Returns
    -------
    actins : (N M,) ndarray of actins_dtype
        List of frames. Each frame is a list of records.

    """
    actins = []
    for i, c in enumerate(coords):
        for x, y in c:
            actins.append((x, y, radius, i))
    return np.array(actins, dtype=actins_dtype)


def motors2coords(motors):
    """
    Convert motor frames to coordinates and attachment locations.

    Parameters
    ----------
    motors : (N,) ndarray of motors_dtype
        List of frames. Each frame is a list of records.

    Returns
    -------
    coords : (N, 2, 2) ndarray of float
        List of motor coordinates.
        Each motor coordinate is a list of head coordinates.
    attached : (N, 2, 2) ndarray of int
        List of motor attachment locations.
        Each motor attachment location
        is a list of head attachment locations.

    """
    coords = []
    attached = []
    for motor in motors:
        x0 = motor["x0"]
        y0 = motor["y0"]
        x1 = x0 + motor["dx"]
        y1 = y0 + motor["dy"]
        coords.append([(x0, y0), (x1, y1)])
        f0 = motor["f_index0"]
        l0 = motor["l_index0"]
        f1 = motor["f_index1"]
        l1 = motor["l_index1"]
        attached.append([(f0, l0), (f1, l1)])
    return np.array(coords), np.array(attached)


def coords2motors(coords, attached=None):
    """
    Convert motor coordinates and attachment locations to frames.

    Parameters
    ----------
    coords : (N, 2, 2) ndarray of float
        List of motor coordinates.
        Each motor coordinate is a list of head coordinates.
    attached : (N, 2, 2) ndarray of int, optional
        List of motor attachment locations.
        Each motor attachment location
        is a list of head attachment locations.

    Returns
    -------
    motors : (N,) ndarray of motors_dtype
        List of frames. Each frame is a list of records.

    """
    if attached is None:
        attached = np.full(np.shape(coords), -1)
    motors = []
    for i, (c, a) in enumerate(zip(coords, attached)):
        (x0, y0), (x1, y1) = c
        (f0, l0), (f1, l1) = a
        motors.append((x0, y0, x1 - x0, y1 - y0, f0, f1, l0, l1))
    return np.array(motors, dtype=motors_dtype)


# -----------------------------------------------------------------------------
# [dtypes]
# numpy dtypes for record lines in AFINES output files

# data/pe.txt
pe_dtype = [
    # time, since there's not header line
    ("time", float),
    # box parameters
    ("xbox", float),
    ("ybox", float),
    ("delrx", float),
    # energies
    ("filament_stretching_energy", float),
    ("filament_bending_energy", float),
    ("filament_exclusion_energy", float),
    ("filament_confinement_energy", float),
    ("motor_stretching_energy", float),
    ("motor_bending_energy", float),
    ("motor_alignment_energy", float),
    ("motor_confinement_energy", float),
    ("xlink_stretching_energy", float),
    ("xlink_bending_energy", float),
    ("xlink_alignment_energy", float),
    ("xlink_confinement_energy", float),
    # virials
    ("filament_stretching_virial_xx", float),
    ("filament_stretching_virial_xy", float),
    ("filament_stretching_virial_yx", float),
    ("filament_stretching_virial_yy", float),
    ("filament_bending_virial_xx", float),
    ("filament_bending_virial_xy", float),
    ("filament_bending_virial_yx", float),
    ("filament_bending_virial_yy", float),
    ("filament_exclusion_virial_xx", float),
    ("filament_exclusion_virial_xy", float),
    ("filament_exclusion_virial_yx", float),
    ("filament_exclusion_virial_yy", float),
    ("filament_confinement_virial_xx", float),
    ("filament_confinement_virial_xy", float),
    ("filament_confinement_virial_yx", float),
    ("filament_confinement_virial_yy", float),
    ("motor_stretching_virial_xx", float),
    ("motor_stretching_virial_xy", float),
    ("motor_stretching_virial_yx", float),
    ("motor_stretching_virial_yy", float),
    ("motor_bending_virial_xx", float),
    ("motor_bending_virial_xy", float),
    ("motor_bending_virial_yx", float),
    ("motor_bending_virial_yy", float),
    ("motor_alignment_virial_xx", float),
    ("motor_alignment_virial_xy", float),
    ("motor_alignment_virial_yx", float),
    ("motor_alignment_virial_yy", float),
    ("motor_confinement_virial_xx", float),
    ("motor_confinement_virial_xy", float),
    ("motor_confinement_virial_yx", float),
    ("motor_confinement_virial_yy", float),
    ("xlink_stretching_virial_xx", float),
    ("xlink_stretching_virial_xy", float),
    ("xlink_stretching_virial_yx", float),
    ("xlink_stretching_virial_yy", float),
    ("xlink_bending_virial_xx", float),
    ("xlink_bending_virial_xy", float),
    ("xlink_bending_virial_yx", float),
    ("xlink_bending_virial_yy", float),
    ("xlink_alignment_virial_xx", float),
    ("xlink_alignment_virial_xy", float),
    ("xlink_alignment_virial_yx", float),
    ("xlink_alignment_virial_yy", float),
    ("xlink_confinement_virial_xx", float),
    ("xlink_confinement_virial_xy", float),
    ("xlink_confinement_virial_yx", float),
    ("xlink_confinement_virial_yy", float),
]

# data/filament_e.txt
filament_e_dtype = [
    ("stretching_energy", float),
    ("bending_energy", float),
    ("f_index", int),
]

# txt_stack/actins.txt
actins_dtype = [
    ("x", float),
    ("y", float),
    ("radius", float),
    ("f_index", int),
]

# txt_stack/links.txt
links_dtype = [
    ("x0", float),
    ("y0", float),
    ("dx", float),
    ("dy", float),
    ("f_index", int),
]

# txt_stack/amotors.txt
# txt_stack/pmotors.txt
motors_dtype = [
    ("x0", float),
    ("y0", float),
    ("dx", float),
    ("dy", float),
    ("f_index0", int),
    ("f_index1", int),
    ("l_index0", int),
    ("l_index1", int),
]
