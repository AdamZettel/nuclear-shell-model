import math
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational

# Precompute factorials divided by 10^n, as in the Fortran COMMON block
FCT = [0.0] * 40
FCT[0] = 1.0
for i in range(1, 40):
    FCT[i] = FCT[i - 1] * (i) / 10.0


def vn02ba(j1, j2, j, m1, m2, m):
    """
    Clebsch-Gordan coefficient <j1,m1,j2,m2 | j,m>.
    Arguments are doubled integers, like in the original Fortran code.
    Example: j1=1 means 1/2, j1=2 means 1, etc.
    """
    # Early rejection checks
    if (m1 + m2) != m:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m) > j:
        return 0.0
    if j > j1 + j2 or j < abs(j1 - j2):
        return 0.0

    # Determine ZMIN, ZMAX
    zmin = 0
    if j < j2 - m1:
        zmin = -j + j2 - m1
    if j < j1 + m2 - zmin:
        zmin = -j + j1 + m2

    zmax = j1 + j2 - j
    if j2 + m2 < zmax:
        zmax = j2 + m2
    if j1 - m1 < zmax:
        zmax = j1 - m1

    cc = 0.0
    for z in range(zmin, zmax + 1, 2):
        ja = z // 2 + 1
        jb = (j1 + j2 - j - z) // 2 + 1
        jc = (j1 - m1 - z) // 2 + 1
        jd = (j2 + m2 - z) // 2 + 1
        je = (j - j2 + m1 + z) // 2 + 1
        jf = (j - j1 - m2 + z) // 2 + 1
        fase = (-1) ** (z // 2)
        f2 = fase

        cc += f2 / (
            FCT[ja - 1]
            * FCT[jb - 1]
            * FCT[jc - 1]
            * FCT[jd - 1]
            * FCT[je - 1]
            * FCT[jf - 1]
        )

    ja = (j1 + j2 - j) // 2 + 1
    jb = (j1 - j2 + j) // 2 + 1
    jc = (-j1 + j2 + j) // 2 + 1
    jd = (j1 + m1) // 2 + 1
    je = (j1 - m1) // 2 + 1
    jf = (j2 + m2) // 2 + 1
    jg = (j2 - m2) // 2 + 1
    jh = (j + m) // 2 + 1
    ji = (j - m) // 2 + 1
    jj = (j1 + j2 + j + 2) // 2 + 1

    f1 = j + 1
    cc = math.sqrt(
        f1
        * FCT[ja - 1]
        * FCT[jb - 1]
        * FCT[jc - 1]
        * FCT[jd - 1]
        * FCT[je - 1]
        * FCT[jf - 1]
        * FCT[jg - 1]
        * FCT[jh - 1]
        * FCT[ji - 1]
        / FCT[jj - 1]
    ) * cc

    return cc / math.sqrt(10.0)


def compare_all(max_j=6, tol=1e-10):
    """
    Compare vn02ba vs sympy clebsch_gordan for all valid (j1,j2,j,m1,m2,m)
    up to max_j (doubled).
    """
    mismatches = []
    for j1 in range(0, max_j + 1):
        for j2 in range(0, max_j + 1):
            for j in range(abs(j1 - j2), j1 + j2 + 1, 2):
                for m1 in range(-j1, j1 + 1, 2):
                    for m2 in range(-j2, j2 + 1, 2):
                        m = m1 + m2
                        if abs(m) > j:
                            continue
                        val_fortran = vn02ba(j1, j2, j, m1, m2, m)
                        j1s, j2s, js = Rational(j1, 2), Rational(j2, 2), Rational(j, 2)
                        m1s, m2s, ms = Rational(m1, 2), Rational(m2, 2), Rational(m, 2)
                        val_sympy = float(clebsch_gordan(j1s, j2s, js, m1s, m2s, ms))
                        if abs(val_fortran - val_sympy) > tol:
                            mismatches.append(((j1, j2, j, m1, m2, m), val_fortran, val_sympy))
    return mismatches


if __name__ == "__main__":
    mismatches = compare_all(max_j=6)
    if mismatches:
        print("Mismatches found:", len(mismatches))
        for entry in mismatches[:10]:
            print(entry)
    else:
        print("All tests passed")
