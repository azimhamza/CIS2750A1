import os; 
import sqlite3
import MolDisplay

# This class is used to store molecules in a database
class Database:
    def __init__(self, reset=False):
        if reset:
            try:
                # Delete the file if it exists
                os.remove('molecules.db')
            except OSError:
                pass
        # Connect to the database
        self.conn = sqlite3.connect('molecules.db')
        # Create the tables if they do not exist
        self.create_tables()

#---------------------------------------------------------------
   
    def create_tables(self):
        # Connect to the database
        c = self.conn.cursor()
        # Create Elements table
        c.execute('''CREATE TABLE IF NOT EXISTS Elements (
                        ELEMENT_NO INTEGER NOT NULL,
                        ELEMENT_CODE VARCHAR(3) NOT NULL PRIMARY KEY,
                        ELEMENT_NAME VARCHAR(32) NOT NULL,
                        COLOUR1 CHAR(6) NOT NULL,
                        COLOUR2 CHAR(6) NOT NULL,
                        COLOUR3 CHAR(6) NOT NULL,
                        RADIUS DECIMAL(3) NOT NULL)''')
        
        # Create Atoms table
        c.execute('''CREATE TABLE IF NOT EXISTS Atoms (
                        ATOM_ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                        ELEMENT_CODE VARCHAR(3) NOT NULL,
                        X DECIMAL(7,4) NOT NULL,
                        Y DECIMAL(7,4) NOT NULL,
                        Z DECIMAL(7,4) NOT NULL,
                        FOREIGN KEY (ELEMENT_CODE) REFERENCES Elements(ELEMENT_CODE))''')
        
        # Create Bonds table
        c.execute('''CREATE TABLE IF NOT EXISTS Bonds (
                        BOND_ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                        A1 INTEGER NOT NULL,
                        A2 INTEGER NOT NULL,
                        EPAIRS INTEGER NOT NULL)''')
        
        # Create Molecules table
        c.execute('''CREATE TABLE IF NOT EXISTS Molecules (
                        MOLECULE_ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                        NAME TEXT NOT NULL UNIQUE)''')
        
        # Create MoleculeAtom table
        c.execute('''CREATE TABLE IF NOT EXISTS MoleculeAtom (
                        MOLECULE_ID INTEGER NOT NULL,
                        ATOM_ID INTEGER NOT NULL,
                        PRIMARY KEY (MOLECULE_ID, ATOM_ID),
                        FOREIGN KEY (MOLECULE_ID) REFERENCES Molecules(MOLECULE_ID),
                        FOREIGN KEY (ATOM_ID) REFERENCES Atoms(ATOM_ID))''')
        
        # Create MoleculeBond table
        c.execute('''CREATE TABLE IF NOT EXISTS MoleculeBond (
                        MOLECULE_ID INTEGER NOT NULL,
                        BOND_ID INTEGER NOT NULL,
                        PRIMARY KEY (MOLECULE_ID, BOND_ID),
                        FOREIGN KEY (MOLECULE_ID) REFERENCES Molecules(MOLECULE_ID),
                        FOREIGN KEY (BOND_ID) REFERENCES Bonds(BOND_ID))''')
        
        # Save changes and close connection
        self.conn.commit()
        

#---------------------------------------------------------------

    def __setitem__(self, table, values):    
        # Connect to the database
        c = self.conn.cursor()
        
        # Set the values in the table based on the key and tuple of values
        if table == 'Elements':
            c.execute('''INSERT OR IGNORE INTO Elements (ELEMENT_NO, ELEMENT_CODE, ELEMENT_NAME, COLOUR1, COLOUR2, COLOUR3, RADIUS)
                            VALUES (?, ?, ?, ?, ?, ?, ?)''', values)
        elif table == 'Atoms':
            c.execute('''INSERT INTO Atoms (ELEMENT_CODE, X, Y, Z)
                            VALUES (?, ?, ?, ?)''', values)
        elif table == 'Bonds':
            c.execute('''INSERT INTO Bonds (A1, A2, EPAIRS)
                            VALUES (?, ?, ?)''', values)
        elif table == 'Molecules':
            c.execute('''INSERT INTO Molecules (NAME)
                            VALUES (?)''', values)
        elif table == 'MoleculeAtom':
            c.execute('''INSERT INTO MoleculeAtom (MOLECULE_ID, ATOM_ID)
                            VALUES (?, ?)''', values)
        elif table == 'MoleculeBond':
            c.execute('''INSERT INTO MoleculeBond (MOLECULE_ID, BOND_ID)
                            VALUES (?, ?)''', values)
        else:
            raise ValueError('Invalid table name')
        
        # Save changes and close connection
        self.conn.commit()
        

#---------------------------------------------------------------

    def add_atom(self, molname, atom):
        # Connect to the database
        c = self.conn.cursor()
        
        # Insert the atom attributes into the Atoms table
        c.execute('''INSERT INTO Atoms (ELEMENT_CODE, X, Y, Z)
                        VALUES (?, ?, ?, ?)''', (atom.element, atom.x, atom.y, atom.z))
        
        # Get the ID of the newly inserted atom
        atom_id = c.lastrowid
        
        # Get the ID of the named molecule
        c.execute("SELECT MOLECULE_ID FROM Molecules WHERE NAME=?", (molname,))
        molecule_id = c.fetchone()[0]
        
        # Insert an entry into the MoleculeAtom table linking the named molecule to the atom entry in the Atoms table
        c.execute('''INSERT INTO MoleculeAtom (MOLECULE_ID, ATOM_ID)
                        VALUES (?, ?)''', (molecule_id, atom_id))
        
        # Save changes and close connection
        self.conn.commit()
        

#---------------------------------------------------------------

    def add_bond(self, molname, bond):
        # Connect to the database
        c = self.conn.cursor()
        
        # Insert the bond attributes into the Bonds table
        c.execute('''INSERT INTO Bonds (A1, A2, EPAIRS)
                        VALUES (?, ?, ?)''', (bond.a1, bond.a2,bond.epairs))
        
        # Get the ID of the newly inserted bond
        bond_id = c.lastrowid
        
        # Get the ID of the named molecule
        c.execute("SELECT MOLECULE_ID FROM Molecules WHERE NAME=?", (molname,))
        molecule_id = c.fetchone()[0]
        
        # Insert an entry into the MoleculeBond table linking the named molecule to the bond entry in the Bonds table
        c.execute('''INSERT INTO MoleculeBond (MOLECULE_ID, BOND_ID)
                        VALUES (?, ?)''', (molecule_id, bond_id))
        
        # Save changes and close connection
        self.conn.commit()

#---------------------------------------------------------------

    def add_molecule(self, name, fp):
    # Import the necessary packages
        
        from MolDisplay import Molecule
        
        # Create a Molecule object and parse its contents
        mol = Molecule()
        mol.parse(fp)
        
        # Connect to the database
        c = self.conn.cursor()
        
        # Add an entry to the Molecules table
        c.execute("INSERT INTO Molecules (NAME) VALUES (?)", (name,))
        
        # Get the ID of the newly inserted molecule
        molecule_id = c.lastrowid
        
        # Add each atom to the Atoms table and link it to the named molecule in the MoleculeAtom table
        i = 0
        while i < mol.atom_no:
            atom = mol.get_atom(i)
            i += 1
            self.add_atom(name, atom)
        
        # Add each bond to the Bonds table and link it to the named molecule in the MoleculeBond table
        j = 0
        while j < mol.bond_no:
            bond = mol.get_bond(j)
            j += 1
            self.add_bond(name, bond)
        
        # Save changes and close connection
        self.conn.commit()

#---------------------------------------------------------------

        
if __name__ == "__main__":
        db = Database(reset=True)
        db.create_tables();

        db['Elements'] = (1, 'H', 'Hydrogen', 'FFFFFF', '050505', '020202', 25)
        db['Elements'] = (6, 'C', 'Carbon', '808080', '010101', '000000', 40)
        db['Elements'] = (7, 'N', 'Nitrogen', '0000FF', '000005', '000002', 40)
        db['Elements'] = (8, 'O', 'Oxygen', 'FF0000', '050000', '020000', 40)

        fp = open('water-3D-structure-CT1000292221.sdf')
        db.add_molecule('Water', fp)
        fp.close()

        fp = open('caffeine-3D-structure-CT1001987571.sdf')
        db.add_molecule('Caffeine', fp)
        fp.close()

        fp = open('CID_31260.sdf')
        db.add_molecule('Isopentanol', fp)
        fp.close()

        # display tables
        print(db.conn.execute("SELECT * FROM Elements;").fetchall())
        print(db.conn.execute("SELECT * FROM Molecules;").fetchall())
        print(db.conn.execute("SELECT * FROM Atoms;").fetchall())
        print(db.conn.execute("SELECT * FROM Bonds;").fetchall())
        print(db.conn.execute("SELECT * FROM MoleculeAtom;").fetchall())
        print(db.conn.execute("SELECT * FROM MoleculeBond;").fetchall())

            
