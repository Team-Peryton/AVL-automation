class Plane():
    def __init__(self,
                name:str,
                geom_file:str=None,
                results_file:str=None,
                Xcg:float=None,
                np:float=None,
                sm:float=None,
                sm_ideal:float=None,
                Lt:float=None,
                St_h:float=None,
                St_v:float=None,
                mac:float=None,
                ARt:float=None,
                span:float=None,
                dihedral_angle:float=None,
                dihedral_split:float=None,
                dihedral_splitY:float=None,
                tipY:float=None,
                tipZ:float=None,
                Xle:float=None,
                cases:list=None,
                polars=None,
                eigen_modes=None,
                tail_config=None,
                b_th=None,
                b_tv=None,
                c_t=None
                ):

        self.name=name
        self.geom_file=geom_file
        self.results_file=results_file
        self.Xcg=Xcg
        self.np=np
        self.sm=sm
        self.sm_ideal=sm_ideal
        self.Lt=Lt
        self.St_h=St_h  #   Equivilent horizontal tail area
        self.St_v=St_v  #   Equivilent vertical tail area
        self.mac=mac
        self.ARt=ARt
        self.span=span
        self.dihedral_angle=dihedral_angle
        self.dihedral_split=dihedral_split
        self.dihedral_splitY=dihedral_splitY
        self.tipY=tipY
        self.tipZ=tipZ
        self.Xle=Xle
        self.cases=cases
        self.polars=polars
        self.eigen_modes=eigen_modes
        self.tail_config=tail_config
        self.b_th=b_th
        self.b_tv=b_tv
        self.c_t=c_t

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
        with open(file,'r') as file:
            for line in file.readlines():
                if line=="\n" or line[0]=="#":  #   Ignores blank lines & comments
                    continue
                else:
                    file_str.append(line)
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
            if line.split()[0]==section_name:
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