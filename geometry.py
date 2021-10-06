class Plane:
    def __init__(self,
                name:str,
                geom_file:str=None,
                results_file:str=None,
                Xcg:float=None,
                np:float=None,
                sm:float=None,
                Lt:float=None,
                St:float=None,
                mac:float=None,
                ARt:float=None,
                d_theta:float=None,
                Xle:float=None,
                cases:list=None,
                polars=None
                ):

        self.name=name
        self.geom_file=geom_file
        self.results_file=results_file
        self.Xcg=Xcg
        self.np=np
        self.sm=sm
        self.Lt=Lt
        self.St=St
        self.mac=mac
        self.ARt=ARt
        self.d_theta=d_theta
        self.Xle=Xle
        self.cases=cases
        self.polars=polars
    
    def make_reference(self,plane_geom:list)->list:
        if self.name=="reference":
            ref_plane=list()
            surface=False
            section=False

            for line in plane_geom:
                if line.split()[0]=="Elevator": #   Finds elevator surface
                    surface=True
                if line.split()[0]=="SECTION" and surface==True:    #   Finds section
                    section=True        
                if section!=True:   #   Adds everything that isn't section 
                    ref_plane.append(line)     
                if line.split()[0]=="SURFACE" and surface==True:    #   Finds next surface
                    ref_plane.append("YES PLEASE\n")    #   Adds marker for adding new sections
                    section=False
                    surface=False
                    ref_plane.append(line)  #   Adds the rest of the file
        else:
            print("Plane must be reference.")
            
        return ref_plane

    def calc_SM(self):
        with open(self.results_file,'r') as text:
            lines=text.readlines()[50]

        self.np=float(lines.split()[-1])
        self.sm=round((self.np-self.Xcg)/self.mac,2)

    def make_dihedral_ref(self,plane_geom:list,split_Yle,aerofoil)->list:
        if self.name=="reference":
            ref_plane_geom=list()
            surface=False
            section=False

            count=0
            for index,line in enumerate(plane_geom):
                if line=="Main Wing\n": #   Finds elevator surface
                    surface=True
                if line.split()[0]=="SECTION" and surface==True:    #   Finds section
                    section=True
                    self.Xle=plane_geom[index+1].split()[0]  #   Gets wing X location  
                    count+=1       
                if count<2 or section==False:   #   Adds everything that isn't section 
                    ref_plane_geom.append(line)     
                if line.split()[0]=="SURFACE" and surface==True:    #   Finds next surface
                    split=Section(self.Xle,split_Yle,0,self.mac,19,-1,aerofoil)  #   Generates split section
                    ref_plane_geom.append("\n"+split.create_input()+"\n")   #   Adds split section
                    ref_plane_geom.append("YES PLEASE\n")    #   Adds marker for adding new sections
                    section=False
                    surface=False
                    ref_plane_geom.append(line)  #   Adds the rest of the file

        return ref_plane_geom

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

    def create_input(self):
        surf_str="{0}\n#Nchord Cspace\n".format(self.name)
        surf_str+="{0} {1}\n".format(self.nchord,self.cspace)

        if self.y_duplicate is not None:
            surf_str+="COMPONENT\n{0}\n".format(self.component)
        surf_str+="YDUPLICATE\n{0}\n".format(self.y_duplicate)
        surf_str+="SCALE\n1 1 1\n"
        surf_str+="TRANSLATE\n0 0 0\n"
        if self.angle is not None:
            surf_str+="ANGLE\n{0}\n".format(self.angle)

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

    def create_input(self):
        section_str="SECTION\n#Xle Yle Zle Chord Ainc Nspan Sspace\n"
        section_str+="{0} {1} {2} {3} {4} {5} {6}\n".format(self.Xle,
                                                            self.Yle,
                                                            self.Zle,
                                                            self.chord,
                                                            0,
                                                            self.nspan,
                                                            self.sspace)
        section_str+="AFIL 0.0.1.0\n{0}\n".format(self.aerofoil)

        return section_str
