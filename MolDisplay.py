import molecule

# -----------------------------Given Code---------------------------------------
radius = 0
element_name = {}
header = """<svg version="1.1" width="1000" height="1000"
xmlns="http://www.w3.org/2000/svg">"""
footer = """</svg>"""
offsetx = 500
offsety = 500
# -------------------------------------------------------------------------------


# Define the Atom class
class Atom:

    # Initialize an Atom object with an atom class/struct as its argument
    def __init__(self, atom):
        # Store the z value of the wrapped atom
        self.atom = atom
        self.z = atom.z
        # Compute the z value based on x and y coordinates

        # Store the wrapped atom itself in the atom attribute

    # Return a string that displays the element, x, y, and z values of the wrapped atom
    def __str__(self):
        return 'element={} x={} y={} z={}'.format(self.atom.element, self.atom.x, self.atom.y, self.z)

    # Return an SVG circle element string representing the atom

    def svg(self):
        # Compute the cx, cy, r, and fill attributes for the circle element
        x = self.atom.x * 100.0 + offsetx
        y = self.atom.y * 100.0 + offsety
        rad = radius[self.atom.element]
        colour = element_name[self.atom.element]

        # Return the SVG circle element string in the desired format
        return '  <circle cx="%.2f" cy="%.2f" r="%d" fill="url(#%s)"/>\n' % (x, y, rad, colour)


# Define the Bond class
class Bond:

    # Initialize a Bond object with a bond class/struct as its argument
    def __init__(self, c_bond):
        self.bond = c_bond
        # Compute the z value based on x and y coordinates of both atoms in the bond
        self.z = c_bond.z
        # Store the wrapped bond itself in the bond attribute

    # Return a string that displays the relevant information of the wrapped bond
    def __str__(self):
        return 'a1.x1={} a1.y1={} a2.x1={} a2.y2={} bond.dx={} bond.dy={} len={} z value={}'.format(
            self.bond.x1, self.bond.y1, self.bond.x2, self.bond.y2, self.bond.dx, self.bond.dy, self.bond.len, self.z)

    # Return an SVG polygon element string representing the bond

    def svg(self):
        # Calculate the coordinates of the 4 corners of the thick line between the atoms
        xa = (self.bond.x1 * 100.0) + offsetx - (self.bond.dy * 10.0)
        ya = (self.bond.y1 * 100.0) + offsety + (self.bond.dx * 10.0)
        xb = (self.bond.x1 * 100.0) + offsetx + (self.bond.dy * 10.0)
        yb = (self.bond.y1 * 100.0) + offsety - (self.bond.dx * 10.0)
        xc = (self.bond.x2 * 100.0) + offsetx + (self.bond.dy * 10.0)
        yc = (self.bond.y2 * 100.0) + offsety - (self.bond.dx * 10.0)
        xd = (self.bond.x2 * 100.0) + offsetx - (self.bond.dy * 10.0)
        yd = (self.bond.y2 * 100.0) + offsety + (self.bond.dx * 10.0)

        # Return the SVG polygon element string in the desired format
        return '  <polygon points="%.2f,%.2f %.2f,%.2f %.2f,%.2f %.2f,%.2f" fill="green"/>\n' % (
            xa, ya, xb, yb, xc, yc, xd, yd)


def JoinString(s):
    return ' '.join(s.split())


class Molecule(molecule.molecule):
    def __str__(self):
        atom_str = "\n".join([str(Atom(self.get_atom(i))) for i in range(self.atom_no)])
        bond_str = "\n".join([str(Bond(self.get_bond(i))) for i in range(self.bond_no)])
        return f"{atom_str}\n{bond_str}"

    def svg(self):
        string = header
        list = []
        i, j = 0, 0

        while i < self.atom_no and j < self.bond_no:
            if self.get_atom(i).z < self.get_bond(j).z:
                list.append(Atom(self.get_atom(i)))
                A = Atom(self.get_atom(i))
                append = A.svg()
                string += append
                i += 1
            else:
                list.append(Bond(self.get_bond(j)))
                B = Bond(self.get_bond(j))
                append = B.svg()
                string += append
                j += 1

        while i < self.atom_no:
            list.append(Atom(self.get_atom(i)))
            A = Atom(self.get_atom(i))
            append = A.svg()
            string += append
            i += 1

        while j < self.bond_no:
            list.append(Bond(self.get_bond(j)))
            B = Bond(self.get_bond(j))
            append = B.svg()
            string += append
            j += 1

        string += footer

        return string

    def parse(self, f):
        for i in range(3):  # skip first 3 lines
            f.readline()

        first_line = JoinString(f.readline().strip())
        first_line_array = first_line.split(' ')
        atom_m = int(first_line_array[0])
        bond_m = int(first_line_array[1])

        for i in range(atom_m):  # parse through atoms
            atom_line = JoinString(f.readline().strip())
            atom_line_array = atom_line.split(' ')
            self.append_atom(atom_line_array[3], float(atom_line_array[0]), float(
                atom_line_array[1]), float(atom_line_array[2]))

        for i in range(bond_m):  # parse through bonds
            bond_line = JoinString(f.readline().strip())
            bond_line_array = bond_line.split(' ')
            # Moved the -1 operation here
            self.append_bond(int(bond_line_array[0]) - 1, int(
                bond_line_array[1]) - 1, int(bond_line_array[2]))

    # Sort the molecule after parsing is complete
    def do_POST(self):
        self.atoms.sort(key=lambda x: x.atom_id)
        self.bonds.sort(key=lambda x: x.bond_id)

   
