class Plane():
    def __init__(self,name:str,geom_file=None):

        self.name=name
        self.geom_file=geom_file
        self.results_file=None
        self.Xcg=None
        self.np=None
        self.sm=None
        self.sm_ideal=None
        self.Lt=None
        self.St_h=None  #   Equivilent horizontal tail area
        self.St_v=None  #   Equivilent vertical tail area
        self.mac=None
        self.ARt=None
        self.span=None
        self.dihedral_angle=None
        self.dihedral_split=None
        self.dihedral_splitY=None
        self.tipY=None
        self.tipZ=None
        self.Xle=None
        self.cases=None
        self.polars=None
        self.modes=None
        self.tail_config=None
        self.b_th=None
        self.b_tv=None
        self.c_t=None
        self.theta=None

        if geom_file!=None:
            self.read(geom_file)

        return None

    def read(self,file:str):
        """
        Loads AVL plane file.

        Parameters:
        -----------
        file: str; AVL plane geometry file.

        Returns:
        -------
        file_str: list; A list of lines in the plane file.
        """
        file_str=list()
        with open(file,'r') as f:
            lines=f.readlines()

        wing=False

        for i,line in enumerate(lines):
            if line=="\n" or line[0]=="#":  #   Ignores blank lines & comments
                continue
            else:
                file_str.append(line)

            if line.strip()=="Main Wing":
                wing=True
            if line.split()[0]=="SECTION" and wing==True:
                Xle=lines[i+1].split()[0]
                wing=False

        self.Xle=Xle

        ref_dims=[n for n in file_str][3]
        self.mac=float(ref_dims.split()[1])
        self.span=float(ref_dims.split()[2])
        self.ARw=self.span/self.mac

        self.file_str=file_str

        return None

    def strip_section(self,section_name):   
        """
        Removes section by name.

        Parameters:
        ----------
        section: str; Name of section in plane file to remove.

        Returns:
        --------
        None - section is stripped inplace.
        """
        stripped_str=[]

        surface=False
        section=False
        found=False
        for line in self.file_str:
            if line.strip()==section_name:
                found=True
                surface=True
            if line.split()[0]=="SECTION" and surface==True:
                section=True
            if section!=True:   #   Adds everything that isn't section 
                stripped_str.append(line)
            if line.split()[0]=="SURFACE" and surface==True:    #   Finds next surface
                stripped_str.append("MARKER\n")    #   Adds marker for adding new sections
                section=False
                surface=False
                stripped_str.append(line)  #   Adds the rest of the file
        
        self.file_str=stripped_str

        if found==False:
            raise KeyError(f"Section '{section_name}' not found.")

        return None

    def strip_surface(self,surface_name):   
        """
        Removes surface by name.

        Parameters:
        ----------
        section: str; Name of surface in plane file to remove.

        Returns:
        --------
        None - surface is stripped inplace.
        """
        stripped_str=[]

        surface=False
        found=False
        for i,line in enumerate(self.file_str):
            if line.split()[0]=="SURFACE":
                if self.file_str[i+1].split()[0]==surface_name:
                    found=True
                    surface=True
                if self.file_str[i+1].split()[0]!=surface_name:
                    surface=False

            if surface==False:
                stripped_str.append(line)

        self.file_str=stripped_str

        if found==False:
            raise KeyError(f"Surface '{surface_name}' not found.")

        return None
    
    def calc_SM(self):
        """
        Reads stability analysis results file and calculates SM based
        on MAC, neutral point, and Xcg.

        Returns:
        sm: float; Static margin.
        """
        with open(self.results_file,'r') as text:
            lines=text.readlines()[50]

        self.np=round(float(lines.split()[-1]),1)
        self.sm=round((self.np-self.Xcg)/self.mac,2)

        return self.sm

    def calc_Xcg_ideal(self):
        """
        Reads stability analysis results file and calculates ideal Xcg
        based on MAC, neutral point, and ideal SM.

        Returns:
        Xcg: float; Ideal CG location in X.
        """
        with open(self.results_file,'r') as text:
            lines=text.readlines()[50]

        self.np=round(float(lines.split()[-1]),1)
        self.Xcg=round(self.np-(self.mac*self.sm_ideal),1)

        return self.Xcg

class Surface():
    #Creates surface (eg wing type)
    def __init__(self,name,nchord,cspace,component,aerofoil,y_duplicate=None,angle=None):
        self.name=name
        self.nchord=nchord
        self.cspace=cspace
        self.component=component
        self.y_duplicate=y_duplicate
        self.angle=angle
        self.aerofoil=aerofoil

    def string(self):
        surf_str    =   f"{self.name}\n#Nchord Cspace\n"
        surf_str    +=  f"{self.nchord} {self.cspace}\n"

        if self.y_duplicate is not None:
            surf_str+=  f"COMPONENT\n{self.component}\n"

        surf_str    +=  f"YDUPLICATE\n{self.y_duplicate}\n"
        surf_str    +=  f"SCALE\n1 1 1\n"
        surf_str    +=  f"TRANSLATE\n0 0 0\n"

        if self.angle is not None:
            surf_str+=  f"ANGLE\n{self.angle}\n"

        return surf_str

class Section():
    def __init__(self,Xle,Yle,Zle,chord,nspan,sspace,aerofoil):
        self.Xle=Xle
        self.Yle=Yle
        self.Zle=Zle
        self.chord=chord
        self.nspan=nspan
        self.sspace=sspace
        self.aerofoil=aerofoil

    def string(self):
        section_str="SECTION\n#Xle Yle Zle Chord Ainc Nspan Sspace\n"
        section_str+=f"{self.Xle} {self.Yle} {self.Zle} {self.chord} {0} {self.nspan} {self.sspace}\n"
        section_str+=f"AFIL 0.0.1.0\n{self.aerofoil}\n"

        return section_str

if __name__=="__main__":
    plane=Plane('aria3','Aria3.avl')
