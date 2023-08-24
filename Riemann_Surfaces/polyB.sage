print("ver 2020.07.24")
#2019.4.25
#Agregue AccionAut que es para encontrar representantes de vectores generadores solo con Aut(G) y no Braid.
#Implemented an output of intermediate_covers and introduced a verbose parameter to suppress printing
#2018.3.16
#Added rational characters involved in group action to Poly
#Addded central idempotents as a list corresponding to each rational character
#2017.3.21
#Reimplemented words as belonging to free group instead of string.
#July 29th 2016
#April 5th 2016
#print "ver 2016.04.05"
#implemented permutation_representation
#fixed bug in split_spaces which caused non symplectic representations in Poly for some groups.
#August 25th 2015
#Improvements of Orbit, AutOrbit, BraidOrbit
#July 25th 2014
#Added as_word to translate from permutation group to finitely presented group.  Needs testing.
#Changed output of riemann_hurwitz to omit empty lists.
#20.12.2013
#Corrected translation error in Orbits function
#Changed find_all_generators function to not sort generating vectors
#Added Orbit, BraidOrbit and AutOrbit.  These still need comments.
#18.10.2013
#added SmallGroup wrapper for gap function
#25.6.2013
#added a change of basis attribute to Poly
#23-11-2012
#Changed edge-labeling algorithm
#8-8-2012
#Reduced the number of equations used to find moebius_invariant so it would run for G=S_5
#Now instead of using all the group elements, it only uses the generators.
#14-5-2012
#Corrected a bug which meant that often the returned symplectic generators were wrong
#In addition, now P.moebius_invariant() should work.

#21-11-2011
#P.set_color is now implemented (P a Poly)
#some improvements made to shorten_words routine
#seed option for Poly so one can choose some of the faces that will be used as class representatives.

#version 17-11-2011
#added routine to shorten words 
#version 10-11-2011
#added g_action2 which returns the image of the loop without moving to the border

#version 29-9-2011

#version 27-9-2011
#added some documentation and fixed problem with non-ascii characters

#version 7-4-2011

#TODO:  1) Improve the better_basis algorithm so it may also work when the group is small relative to the genus.
#       2) Work on getting a convex polygon. Can it be done always?
#       3) Reduce words such as aa,b^-1b^-1,a^0, etc DONE 
#       4) Allow flexibility in the construction of the fundamental polygon covering the sphere
#       5) Add verifications 
i=CDF(I) #define i as the complex number i
Pi=RDF(pi) #define the variable Pi as the real value pi

libgap.eval('LoadPackage("kbmag")')

class Vertex():
    def __init__(self,pw=0,fuchsian_element=None,xword=None,base_stabilizer=0):
        self.label=None
        self.weight=1
        self.potential_weight=pw
        self.xw=xword
        self.base=base_stabilizer
        self.w=fuchsian_element  #stabilizer
    def __repr__(self):
        return 'V%d'%self.label
    def _latex_(self):
        return 'V%d'%self.label
    def is_convex(self):
        return False if 2*self.weight>self.potential_weight and self.weight<self.potential_weight else True

class GeodesicArc(list):
    def __init__(self, begin,end,rgbcolor=[0,0,1],**keywords):
        """
        Given two complex numbers inside the unit disc, this defines the hyperbolic arc between them.  The color is blue by default.
        Since in its euclidean representation, hyperbolic geodesics are circles, it is useful to know the center of this euclidean circle.
        It can be obtained as A.center() when A is a geodesic arc.

        INPUT:
                ``begin`` - a complex number z with |z|<1.
                ``end`` - a complex number w with |w|<1.
                ``rgbcolor`` -         

        EXAMPLE::
                sage: A=GeodesicArc((1+i)/3,(1-2*i)/4)
                sage: A.draw()
                sage: A.center()
        """
        list.__init__(self,[begin,end])
        self.color=rgbcolor
        self._center=None
        self.epsilon=0.01
    def A(self):
        """
        Returns the complex value of its initial point on the Poincare disc
        """
        return self[0]
    def B(self):
        """
        Returns the complex value of its ending point on the Poincare disc
        """
        return self[1]
    def midpoint(self):
        """
        Returns the complex value of the midpoint of the geodesic arc.
        """
        return self.isosceles(2)
    def center(self):
        """
        Returns the center of the euclidean circle containing the geodesic arc or 0 if the geodesic arc is too close to being a
        euclidean line.
        """
        if self._center!=None:
            return self._center
        m=matrix(RDF,2,[2*real(self[0]),2*imag(self[0]),2*real(self[1]),2*imag(self[1])])
        x,y=var('__x,__y')
        v=matrix(RDF,2,1,[abs(self[0])^2+1,abs(self[1])^2+1])
        if abs(m.determinant())>self.epsilon:
            S=m.inverse()*v
            self._center=CDF(S[0][0]+i*S[1][0])
            a=arg(self[0]-self._center)
            b=arg(self[1]-self._center)
            if abs(a-b)<0.1:
                self._center=CDF(0)
            return self._center #complex number corresponding to the center of the euclidean circle for the hyperbolic arc
        else:
            self._center=CDF(0) 
            return CDF(0)       
    def __repr__(self):
        return "Segment from "+str(self[0])+ " to "+str(self[1])
    def draw(self,**keywords):
        """
        Return the graphical object consisting of a line through the points approximating this geodesic arc
        """
        C=self.center()
        if C == CDF(0):
            return line([[real(self[0]),imag(self[0])],[real(self[1]),imag(self[1])]],rgbcolor=self.color,**keywords)
        a=arg(self[0]-C)
        b=arg(self[1]-C)
        r=RDF(abs(self[0]-C))
        if a>b:
            c=a
            a=b
            b=c
        if a*b<0:
            if b-a>Pi:
                c=a+2*Pi
                a=b
                b=c
        x=real(C)
        y=imag(C)
        return arc((x,y),r,sector=(a,b),rgbcolor=self.color,**keywords)
    def set_color(self,rgbcolor):
        """
        Sets the color of this geodesic arc.  Expects rgb input
        """
        self.color=rgbcolor
    def rotate(self,alpha,z_0):
        """
        Rotates the geodesic arc around the complex point z_0 in a counterclockwise angle alpha.
        sage: p1=GeodesicArc((1+i)/2,(1-2*i)/5)
        sage: p1.rotate(Pi/2,0)
        sage: p1
        """
        R=rotation_matrix(alpha,z_0) 
        self.apply_moebius(R)
    def rotated(self,alpha,z_0):
        """
        Returns a new geodesic arc as a result of rotating the current one around the complex point z_0 in counterclockwise angle alpha.
        sage: p1=GeodesicArc((1+i)/2,(1-2*i)/5)
        sage: p2=p1.rotated(Pi/2,0)
        sage: p1,p2
        """
        c=deepcopy(self)    
        c.rotate(alpha,z_0)
        return c
    def move(self,z_0):
        """
        Translates the current one along the vector from 0 to z_0.
        sage: p1=GeodesicArc((1+i)/2,(1-2*i)/5)
        sage: p1.move(0.3-0.2*i)
        sage: p1
        """
        R=matrix(2,[1,z_0,z_0.conjugate(),1])
        self.apply_moebius(R)
    def moved(self,z_0):
        """
        Returns a new geodesic arc as a result of translating the current one along the vector from 0 to z_0.
        sage: p1=GeodesicArc((1+i)/2,(1-2*i)/5)
        sage: p2=p1.moved(0.3-0.2*i)
        sage: p1,p2
        """
        c=deepcopy(self)    
        c.move(z_0)
        return c
    def apply_moebius(self,R):
        self[0]=moebius_transformation(self[0],R)
        self[1]=moebius_transformation(self[1],R)
        self._center=None
    def applied_moebius(self,R):
        c=deepcopy(self)    
        c.apply_moebius(R)
        return c
    def reflect(self):
        self[0]=self[0].conjugate()
        self[1]=self[1].conjugate()
    def reflected(self):
        c=deepcopy(self)
        c.reflect()
        return c
    def length(self):
        """
        Returns the hyperbolic length of the arc.
        
        Example:
        sage: A=GeodesicArc((1+i)/2,(1-2*i)/5)
        sage: A.length()
        """
        return 2*asinh(sqrt((self[0]-self[1]).abs()^2/((1-self[0].abs()^2)*(1-self[1].abs()^2))))
    def isosceles(self,n):
        alpha=2*Pi/n
        b=asinh(sqrt((cosh(self.length())-1)/(1-cos(alpha))))
        beta=asin(sinh(b)*sin(alpha)/sinh(self.length()))
        z=(self[0]-self[1])/(1-self[0]*self[1].conjugate())
        w=z/z.abs()*tanh(b/2)*exp(i*beta)
        return (w+self[1])/(1+w*self[1].conjugate())
    def __add__(self,x):
        return GeodesicFigure([self])+GeodesicFigure([x])
    def __radd__(self,x):
        if x==0:
            return self
    def __cmp__(self,other):
        if list(self)==list(other):
            return 0

class Edge(GeodesicArc):
    def __init__(self,face,k,vertices):
        self.rev=[]
        self.face=face
        self.direction=k
        self.label=None
        self.check=False
        GeodesicArc.__init__(self,face[k],face[(k+1)%len(face)])
        self.is_border=None
    def a(self):
        return self.face.vertices[self.direction-1]
    def b(self):    

        return self.face.vertices[self.direction]
    def _left_rotate(self):
        return self.face._edges[self.direction-1].rev
             
    def __repr__(self):
        l=self.label
        s='E%d'%l if l>0 else '-E%d'%-l
        return s

        
    def __str__(self):
        return str(self.label)
    
    def _latex_(self):
        l=self.label
        s='E%d'%l if l>0 else '-E%d'%-l
        return s
        
    def costarters(self):
        next=self._left_rotate()
        r=[self]
        while next!=self:
            r.append(next)
            next=next._left_rotate()
        return r
    def coterminals(self):
        return[k.rev for k in self.rev.costarters()]
        
    def draw(self,label=False,fontsize=10,**keywords):
        l=GeodesicArc.draw(self,**keywords)
        if label=='edges':
            C=(self[0]+self[1])/2
            l=l+text(str(self.label),(C.real(),C.imag()),fontsize=fontsize)
        if label=='vertices':
            C=self[0]
            l=l+text(str(self.a().label),(C.real(),C.imag()),fontsize=fontsize)
        return l

class Path(list):
    def __init__(self, edges):
        """
        Returns a path consisting of the given edges in the given order.  If the edges do not connect it produces an error.
        """
        p=edges[0].b()
        for k in edges[1:]:
            if not type(k)==type(edges[0]):
                raise ValueError("not all edges are of the same type")
            if not k.a().label==p.label:
                raise ValueError("edges do not connect")
            p=k.b()
        list__init__(self, edges)

class GeodesicPolygon(list):
    def __init__(self, points,label='',rgbcolor=[0,0,1]):
        """
        Returns the hyperbolic polygon with vertices at the points given by complex numbers inside the unit disc
        
        EXAMPLE::
        
            sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
            sage: p1.draw()
        
        INPUT::
        
            -``points`` - a list of complex numbers of modulus less than 1.
            -``label`` - a string (optional)
            -``rgbcolor`` - rgbvalue for the color (default is [0,0,1] which is blue)
        
        EXAMPLES::
        
            sage: p2=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5],rgbcolor=[1,0,0])
            sage: p2.draw()
            
        Or you could change the color of an existing polygon:
        
        EXAMPLE::
        
            sage: p2.set_color([0,1,0])
            sage: p2.draw()
            
        Other implemented actions are: move, rotate, apply_moebius, reflect, add_isosceles.
        """ 
        for k in points:
            if not k in i.parent():
                raise ValueError("vertices should all be complex numbers")
        if max([k.abs() for k in points])>1:
            raise ValueError("Some vertex is not inside the Poincare disc") 
        list.__init__(self, points)
        self._edges=[GeodesicArc(self[k],self[(k+1)%len(points)],rgbcolor) for k in range(len(self))]
        self._center=sum(points)/len(points)
        self.label=label
    def center(self):
        return self._center
    def __repr__(self): 
        return "Closed Polygon with "+str(len(self))+" vertices"
    def __add__(self,x):
        return GeodesicFigure(self._edges)+GeodesicFigure(x._edges)
    def __radd__(self,x):
        if x==0:
            return self
    def draw(self,label=False,**keywords):
        """
        Returns a graphical object representing the polygon.
        If the polygon's labels are to be displayed, label should be set to true.
        """
        p=sum([self._edges[k].draw(**keywords) for k in range(len(self._edges))])
        p.set_aspect_ratio(1)
        if label:
            p=p+text(self.label,(self.center().real(),self.center().imag()))
        return p
    def set_color(self,rgbcolor):
        """
        Resets the color of all the edges of the polygon.
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p1.set_color([1,0,0])
        sage: p1.draw()
        An individual edge can be recoloured also:
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p1.edges()[0].set_color([1,0,0])
        sage: p1.draw()
        """
        for k in self._edges:
            k.set_color(rgbcolor)
    def rotate(self,alpha,z_0):
        """
        Rotates the polygon around the complex point z_0 in a counterclockwise angle alpha.
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p1.rotate(Pi/3,0.2)
        sage: p1.draw()
        """
        R=rotation_matrix(alpha,z_0)
        self.apply_moebius(R)
    def rotated(self,alpha,z_0):
        """
        Returns a new polygon as a result of rotating the current one around the complex point z_0 in counterclockwise angle alpha.
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p2=p1.rotated(Pi/3,0.2)
        sage: p2.draw()
        """
        c=deepcopy(self)    
        c.rotate(alpha,z_0)
        return c        
    def move(self,z_0):
        """
        Moves the polygon "adding" z_0 (a complex number inside the unit disc) in a hyperbolic manner.
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p1.move(0.2+0.3*i)
        sage: p1.draw()
        """
        R=matrix(2,[1,z_0,z_0.conjugate(),1])
        self.apply_moebius(R)
    def moved(self,z_0):
        """
        Returns a new polygon as a result of moving the current one.
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p2=p1.moved(0.2+0.3*i)
        sage: p1.draw()+p2.draw()
        """
        c=deepcopy(self)    
        c.move(z_0)
        return c
    def apply_moebius(self,R,label=''):
        for k in range(len(self)):
            self[k]=moebius_transformation(self[k],R)
            self._edges[k].apply_moebius(R)
        self._center=sum(self)/len(self)
    def applied_moebius(self,R):
        c=deepcopy(self)    
        c.apply_moebius(R)
        return c
    def reflect(self):
        for k in range(len(self)):
            self[k]=self[k].conjugate()
            self._edges[k].reflect()
    def reflected(self):
        c=deepcopy(self)
        c.reflect()
        return c
    def add_isosceles(self,L):
        Q=GeodesicPolygon([self._edges[k].isosceles(L[k]) for k in range(len(self._edges))])
        R=[]
        for k in range(len(self._edges)):    
            R.append(self[k])
            R.append(Q[k])
        return GeodesicPolygon(R)

class FundamentalPolygon(GeodesicPolygon):
    def __init__(self,A0,label='',beta=0,alphabet="abcdefghijklmnopqrstuvwxyz"):
        self.signature=A0
        n=len(A0)-1
        A=A0[:n]
        P0=CDF(tanh(rho(A0)/2))
        P=regular_polygon(n,P0)
        R1=P.add_isosceles(A)
        R2=R1.moved(-R1[1])    
        alpha=-R2[2].arg()+beta
        R3=R2.rotated(alpha,0*i)
        GeodesicPolygon.__init__(self,R3,label)
        self.rotation_matrices=[]
        for k in range(n):
            self.rotation_matrices.append(rotation_matrix(2*Pi/A[k],self[2*k+1]))
        self.rotation_matrices.append(prod(self.rotation_matrices)^-1)
        self._free_group=Groups().free(names=list(alphabet[:n+1]))
        self._word_reducer=None

    def free_group(self):
        return self._free_group
    def word_reducer(self):
        if self._word_reducer==None:
            self._word_reducer=kbword_reducer(self)
        return self._word_reducer

    def letters(self):
        return self.free_group().gens()
   
    def apply_word_to(self,word,obj):
        """
        A word is interpreted as a Moebius transformation and this transformation is applied to the object given.
        The transformation is made in place, so nothing is returned and the object is changed.
        
        INPUT:
            
            ``word`` - a string composed of letters in self.letters().
            
            ``object`` - a GeodesicArc, GeodesicPolygon or GeodesicFigure
        
        EXAMPLES::
            
            sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
            sage: self.free_group().inject_variables()
            sage: self.apply_word_to(a,p1)
                
        """
        R=word(self.rotation_matrices)
        obj.apply_moebius(R)
    def select_by_words(self,L,S=None):
        """
        Each word in L is interpreted as a Moebius transformation and applied to S.  The resulting copies of S are combined and returned as a GeodesicFigure.  If S is not given, self will be used instead.
        
        INPUT::
            
            ``L`` - a list of words in the free group self.free_group()
            ``S`` - an optional GeodesicFigure to which each of the transformations wil be applied. 
                    If none is supplied, the FundamentalPolygon itself will be used.
        
        EXAMPLES::
            
            sage: P=FundamentalPolygon([4,4,3])
            sage: P.free_group().inject_variables()
            sage: C=P.select_by_words([a,a*b^-1])
            sage: C.draw()
     
        ::
            
            sage: P=FundamentalPolygon([4,4,3])
            sage: p1=GeodesicPolygon([0,0.4,0.08+0.23*i,0.2-0.14*i,0.32+0.23*i])
            sage: C=P.select_by_words([a,a*b^-1],p1) 
                
        """
        if S==None: S=self
        picture=GeodesicFigure([deepcopy(S) for k in L])
        for k in range(len(L)):
            self.apply_word_to(L[k],picture[k])
        return picture

class Face(GeodesicPolygon):
    def __init__(self,fund_poly,g,label,w):
        A=fund_poly.signature
        self.signature=A
        self.g=g
        self.label=label
        reducer=fund_poly.word_reducer()
        self.word=reducer(w)
        m=fund_poly.rotation_matrices[0].parent()(w(fund_poly.rotation_matrices))
        p=[moebius_transformation(k,m) for k in fund_poly]
        GeodesicPolygon.__init__(self,p,label)
        self.vertices=[]
        for k in range(len(fund_poly)):
            pwk=A[-1]*(len(A)-1) if k%2 else A[k//2]
            f_ek=None if k%2 else  w*fund_poly.free_group().gens()[k//2]*w.inverse()
            bsk=w.parent().gens()[len(A)-1 if k%2 else k//2]
            self.vertices.append(Vertex(pw=pwk,fuchsian_element=f_ek,xword=w,base_stabilizer=bsk))
        self._edges=[Edge(self,k,self.vertices) for k in range(len(fund_poly))]
    def __repr__(self):
        return 'F%d'%self.label

           
class GeodesicFigure(list):
    def __init__(self, segments,label=''):
        """
        As a more general construction than a polygon, a GeodesicFigure may be constructed from a list of geodesic arcs (segments) or adding existing polygons to each other.
        sage: p=[(0,0), (0.02967,0.0205), (0.05131,0.01541), (0.06467,0.00905), (0.09457,0.01541), (0.10538,0.01732), (0.10538,0.02177), (0.09711,0.03068), (0.07803,0.0364), (0.07421,0.04658), (0.08057,0.06249), (0.09457,0.07521), (0.11302,0.08475), (0.11174,0.10384), (0.12701,0.11275), (0.15628,0.10002), (0.18427,0.07521), (0.23708,0.03577), (0.24853,0.02686), (0.25744,0.02431), (0.25298,0.04276), (0.25425,0.06885), (0.26125,0.09048), (0.27716,0.11656), (0.2937,0.13056), (0.31278,0.14201), (0.3236,0.14519), (0.3376,0.17446), (0.35477,0.19227), (0.37195,0.20881), (0.38786,0.21963), (0.40885,0.23235), (0.42857,0.23871), (0.448287736084 + 0.258819045103*i), (0.44257,0.21581), (0.43939,0.18909), (0.43812,0.15474), (0.43939,0.12611), (0.4413,0.11593), (0.4511,0.10212), (0.46361,0.08058), (0.46951,0.06009), (0.46917,0.03577), (0.45284,0.03299), (0.43721,0.03542), (0.42054,0.0448), (0.40247,0.05592), (0.38684,0.05904), (0.36947,0.05696), (0.35558,0.05245), (0.33647,0.03716), (0.35002,0.02917), (0.36392,0.02535), (0.3792,0.02292), (0.39969,0.02257), (0.41255,0.02014), (0.42331,0.01562), (0.42957,0.00902), (0.43443,0.00208), (0.43617,0)]
        sage: p=[CDF(k) for k in p]
        sage: p=[GeodesicArc(p[k],p[k+1]) for k in range(len(p)-1)]
        sage: p=GeodesicFigure(p)
        sage: Q=tesselate_disc([3,4,4],klein=True,S=p,bound=0.03,alpha=-Pi/6)
        sage: A=Q.draw()
        sage: A.show(axes=False)

        """
        list.__init__(self, segments)
    def __repr__(self): 
        s='['+', '.join([latex(k)for k in self])+']'
        return s
    def __add__(self,x):
        return GeodesicFigure(list.__add__(self,x))
    def __mul__(self,n):
        return GeodesicFigure(list.__mul__(self,n))
    def draw(self,**keywords):
        p=sum([k.draw(**keywords) for k in self])
        p.set_aspect_ratio(1)
        return p
    def set_color(self,color):
        for k in self:
            k.set_color(color)
    def rotate(self,alpha,z_0):
        """
        Rotates the figure around the complex point z_0 in a counterclockwise angle alpha.
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p2=GeodesicPolygon([1/2,(1+i)/3,(1+2*i)/5])
        sage: f=p1+p2
        sage: g=f.draw()
        sage: f.rotate(Pi/2,0)
        sage: f.set_color([1,0,0])
        sage: g+f.draw()
        """
        R=rotation_matrix(alpha,z_0)
        self.apply_moebius(R)
    def rotated(self,alpha,z_0):
        """
        Returns a new figure as a result of rotating the current one around the complex point z_0 in counterclockwise angle alpha.
        sage: p1=GeodesicPolygon([0,(1+i)/2,(1-2*i)/5])
        sage: p2=GeodesicPolygon([1/2,(1+i)/3,(1+2*i)/5])
        sage: f=p1+p2
        sage: g=f.rotated(Pi/2,0)
        sage: g.set_color([1,0,0])
        sage: (g+f).draw()
        """
        c=deepcopy(self)    
        c.rotate(alpha,z_0)
        return c        
    def move(self,z_0):
        R=matrix(2,[1,z_0,z_0.conjugate(),1])
        self.apply_moebius(R)
    def moved(self,z_0):
        c=deepcopy(self)    
        c.move(z_0)
        return c
    def apply_moebius(self,R,label=''):
        for k in range(len(self)):
            self[k].apply_moebius(R)
    def applied_moebius(self,R):
        c=deepcopy(self)    
        c.apply_moebius(R)
        return c
    def reflect(self):
        for k in range(len(self)):
            self[k].reflect()
    def reflected(self):
        c=deepcopy(self)
        c.reflect()
        return c            
    def _latex_(self):
        s='['+', '.join([latex(k)for k in self])+']'
        return s

class GeodesicRegion(GeodesicFigure):
    def __init__(self,faces,generators):
        GeodesicFigure.__init__(self,[])
        A=[k.order() for k in generators]
        self.fundamental_polygon=FundamentalPolygon(A)
        self.generators=generators
        self._expanded_generators=sum([[k,k^-1]for k in self.generators],[])
        self.face_id={}
        L=self.fundamental_polygon
        self._ewg=ewg=sum([[k,k^-1]for k in L.letters()[:-1]],[])
        self._addface(faces)
        self._gfaces=[k.g for k in self]
        self.__make_border()
    def _addface(self,list_of_new_faces,unsafe=[]):
        eg=self._expanded_generators[:-2]
        L=self.fundamental_polygon
        for nf in list_of_new_faces:
            for k in nf._edges:
                k.is_border=True
            self.face_id[nf.g]=len(self)
            self.append(nf)
            for k in range(len(eg)):
                h=k+(-1 if k%2 else 1)
                cn=nf.g*eg[h]
                w=nf.word*self._ewg[h]
                try:
                    of=self[self.face_id[cn]]
                except(KeyError):
                    pass
                else:
                    tf=Face(L,cn,len(self),w)
                    of._edges[k].rev=nf._edges[h]
                    nf._edges[h].rev=of._edges[k]
                    if (tf.center()-of.center()).abs()<0.0001:
                        if nf.vertices[h]!=of.vertices[k-1]:
                            nf.vertices[h]=of.vertices[k-1]
                            nf.vertices[h].weight+=1
                        if nf.vertices[h-1]!=of.vertices[k]:
                            nf.vertices[h-1]=of.vertices[k]
                            nf.vertices[h-1].weight+=1
                        of._edges[k].is_border=False
                        nf._edges[h].is_border=False
                    else:
                        of._edges[k].is_border=True
                        nf._edges[h].is_border=True

    def __make_border(self):
        eg=self._expanded_generators[:-2]
        current_face=None
        for F in self:
            for E in F._edges:
                if E.is_border==True:
                    current_face=F
                    break
        ng=len(eg)
        first_face=current_face
        k=ng-1
        while not current_face._edges[k].is_border: 
            k-=1
        while current_face._edges[k].is_border: 
            k-=1
            if k<0:
                break 
        current_edge=first_edge=k+1
        pre_border=[]
        ready=False
        while not ready:
            if current_face._edges[current_edge].is_border:
                pre_border.append(current_face._edges[current_edge])
            current_edge=(current_edge+1)%ng #rotate 1 counterclockwise on current face
            if not current_face._edges[current_edge].is_border:
                current_face=current_face._edges[current_edge].rev.face
                current_edge=current_edge+(-1 if current_edge%2 else 1)
            if current_face==first_face and current_edge==first_edge:
                ready=true
        self._border=GeodesicFigure(pre_border)
    def border(self):
        return self._border

class Poly():
    def __init__(self,G,x,seed=[],alphabet="abcdefghijklmnopqrstuvwxyz",compute=True):
        tt=cputime()
        if type(G)==type(WeylGroup(['A',1])):
            G=PermutationGroup(gap.IsomorphismPermGroup(G).Image().GeneratorsOfGroup())
        self._G=G
        if type(x[0])==type(3):
            self.signature=x
            self.generators=find_generating_vector(self._G,x)
        if x[0].parent()==self._G:
            self.generators=list(x)
            if prod(self.generators).order()!=1:
                self.generators.append(prod(self.generators)^-1)
            self.signature=[g.order() for g in self.generators]
        if self._G.subgroup(self.generators).order()!=self._G.order():
            raise ValueError("The given generators do not generate the group")         
        self.fundamental_polygon=FundamentalPolygon(self.signature,alphabet=alphabet)
        self.face_id={}
        self._edges=[]
        self.__process_group(seed)
        self._cs={} #dictionary of costarting edges to store information and save calculation time
        self._im=0  #will eventually contain the intersection matrix
        self._homology=0
        self._sgl=0
        self._sgg=0
        self._rel=0
        self.__make_border()
        self.__make_numbering()
        self._border[0].set_color([1,0,0])
        self.intersection_matrix()
        self._central_idempotents=None
        self._rational_characters=None
        self._pp=None
        if compute:
                self.__compute_symplectic_group()
        #print 'calculations took %f seconds of cputime'%cputime(tt)

    def QG(self):
        return GroupAlgebra(self._G,QQ)   
    def hom(self):
        if self._pp==None:
            Ms=self.symplectic_group_generators()[0].parent()
            h1=dict()
            for f in self.faces:
                h1[f.g]=Ms(f.word(self.symplectic_group_generators()))
            self._pp=self.QG().module_morphism(on_basis=lambda x: h1[x],codomain=Ms)
        return self._pp

    def letters(self):
        return self.fundamental_polygon.letters()
    def border(self):
        return self._border

    def homology(self):
        return self._homology
    
    def free_group(self):
        return self.fundamental_polygon.free_group()
        
    def genus(self):
        return self._Genus

    def __repr__(self):
        return "Hyperbolic polygon for " + str(self._G) + " acting on a curve of genus "+str(self._Genus)+" with signature "+str(self.signature)+" and specific generating vector "+str(self.generators)

    def _latex_(self):
        return "Hyperbolic polygon for " + str(self._G) + " acting on a curve of genus "+str(self._Genus)+" with signature "+str(self.signature)+" and specific generating vector "+str(self.generators)
    def draw(self,S=0,label=False,fontsize=10,**keywords):
        if S==0: S=self._G
        if label == 'faces':
            return sum(f.draw(label=True) for f in self.faces)
        if label:
            return GeodesicFigure([k for k in self.border()]).draw(label=label,fontsize=fontsize,**keywords)+GeodesicFigure([k for k in self._inner_edges]).draw(label=False,**keywords)
        return GeodesicFigure([k for k in self._edges if True or k.face.g in S]).draw(**keywords)
    
    def set_color(self,color):
        for k in self.faces:
            k.set_color(color)
    def edges(self):
        return self._edges
    def group(self):
        return self._G
    def mgroup(self):
        return self.fundamental_polygon.rotation_matrices[0].parent()
    def __process_group(self,seed):
        A=self.signature
        L=self.fundamental_polygon
        n=len(A)
        #Creates the list of generators adding their inverses
        #There are several representations: group elements, moebius matrices, fuchsian pre-images
        #no entry is created for the last generator, effectively ignoring it
        self._expanded_generators=eg=sum([[k,k^-1]for k in self.generators[:-1]],[])
        self._ewg=ewg=sum([[k,k^-1]for k in self.letters()[:-1]],[])
        last_element=self._G.identity()
        ng=len(ewg)
        #
        #setup of central part generated by first generator
        g0=self.generators[0]
        w0=ewg[0]
        N=A[0]
        self.faces=faces=[]
        if seed==[]:
            seed=[self.free_group()(1)]
        for w1 in seed:
            g1=self.group()(w1(self.generators))
            for k in range(N):
                nf=Face(L,g0^k*g1,len(faces),w0^k*w1)
                self._addface([nf])
        ready=False
        while not ready:
            ready=True
            for k in range(ng):
                h=k+(-1 if k%2 else 1)
                for g in range(N):
                    if faces[g]._edges[k].is_border==None:
                        ready=False
                        cn=faces[g].g*eg[k]
                        w=faces[g].word*ewg[k]
                        cf=Face(L,cn,len(faces),w)
                        self._addface([cf])
            N=len(faces)
    
    def _addface(self,list_of_new_faces,unsafe=[]):
        eg=self._expanded_generators
        L=self.fundamental_polygon
        for nf in list_of_new_faces:
            self.face_id[nf.g]=len(self.faces)
            self.faces.append(nf)
            for k in range(len(eg)):
                h=k+(-1 if k%2 else 1)
                cn=nf.g*eg[h]
                w=nf.word*self._ewg[h]
                try:
                    of=self.faces[self.face_id[cn]]
                except(KeyError):
                    pass
                else:
                    tf=Face(L,cn,len(self.faces),w)
                    of._edges[k].rev=nf._edges[h]
                    nf._edges[h].rev=of._edges[k]
                    if (tf.center()-of.center()).abs()<0.0001:
                        if nf.vertices[h]!=of.vertices[k-1]:
                            nf.vertices[h]=of.vertices[k-1]
                            nf.vertices[h].weight+=1
                        if nf.vertices[h-1]!=of.vertices[k]:
                            nf.vertices[h-1]=of.vertices[k]
                            nf.vertices[h-1].weight+=1
                        of._edges[k].is_border=False
                        nf._edges[h].is_border=False
                    else:
                        of._edges[k].is_border=True
                        nf._edges[h].is_border=True
        
    def _remove_face(self,nf):
        eg=self._expanded_generators
        L=self.fundamental_polygon
        n=nf.label
        self.face_id.pop(nf.g)
        self.faces.pop(n)
        for k in range(n,len(self.faces)):
            self.faces[k].label-=1
        for k in range(len(eg)):
            h=k+(-1 if k%2 else 1)
            cn=nf.g*eg[h]
            w=nf.word*self._ewg[h]
            nf.vertices[k].weight-=1
            try:
                of=self.faces[self.face_id[cn]]
            except(KeyError):
                pass
            else:
                tf=Face(L,cn,len(self.faces),w)
                of._edges[k].rev=[]
                of._edges[k].is_border=None
        
    def __make_border(self):
        ng=len(self._expanded_generators)
        L=self.fundamental_polygon
        current_face=self.faces[-1] #face which was added last
        k=ng-1
        while not current_face._edges[k].is_border: 
            k-=1
        while current_face._edges[k].is_border: 
            k-=1 
        current_edge=first_edge=k+1
        pre_border=[]
        ready=False
        while not ready:
            if current_face._edges[current_edge].is_border:
                pre_border.append(current_face._edges[current_edge])
            current_edge=(current_edge+1)%ng #rotate 1 counterclockwise on current face
            if not current_face._edges[current_edge].is_border: 
                current_face=current_face._edges[current_edge].rev.face
                current_edge=current_edge+(-1 if current_edge%2 else 1)
            if current_face==self.faces[-1] and current_edge==first_edge:
                ready=true
        self._border=GeodesicFigure(pre_border)
            
    def __make_numbering(self):
        #the next part assigns numbers to the border edges in such a way that n is identified with -n    
        for j in range(len(self.border())/2):
            if not self.border()[2*j].label:
                self._edges.append(self.border()[2*j])
                self.border()[2*j].label=j+1
                self.border()[2*j].rev.label=-j-1
            else:
                print(j)
        #print [k.label for k in self.border()]
        inner=[edg for f in self.faces for edg in f._edges if not edg.is_border]
        self._border_edges=[None]
        self._border_edges.extend(self._edges)
        self._border_edges.extend(reversed(self._edges))
        self._edges.extend([edg.rev for edg in self._edges])
        self.inner=inner
        j+=2
        for edg in inner:
            if not edg.label:
                self._edges.append(edg) 
                edg.label=j
                edg.rev.label=-j
                j+=1
        self._inner_edges=[k for k in inner if k.label>0]
        #the next part labels the vertices
        j=0
        k=0
        for edg in self.border():
            edg.a().position=k
            k+=1        
            if edg.a().label==None:
                for x in edg.costarters(): 
                    x.a().label=j
                j+=1
        self.number_of_vertices=j
        for edg in inner:
            if edg.a().label==None:
                edg.a().label=j
                j+=1    
        
    #nth homology basis element as list of edges in appropriate order
    def _homology_basis_element(self,n):
        BL=len(self.border()) 
        v=self.homology()[n]
        rr=[]
        ro=[]
        for k in range(BL//2):
            ed=self._edges[k]
            if v[k]==1:
                rr.append(ed)
            if v[k]==-1:
                rr.append(ed.rev)
        ro.append(rr.pop())
        j=0
        while rr:
            j=j%len(rr)
            if ro[-1].b().label==rr[j].a().label:
                ro.append(rr.pop(j))
                j=(j-1)
            else:
                j=(j+1)
        return ro
        
    def equivalent_path(self,pth):
        """
        Given a closed path with at least one vertex on the border, 
        this function finds an equivalent path lying completely on the border.
        """
        rr=[]
        k=Zmod(len(pth))(0)
        #searches for an edge in pth whose starting vertex lies on the border
        while not pth[k].a().label<self.number_of_vertices:
            k+=1
        k0=k
        while true:
            if pth[k].is_border:
                rr.append(pth[k])
                k+=1
                if k==k0: break
                continue
            startingpoint=pth[k].a()
            while not pth[k].b().label<self.number_of_vertices:
                k+=1
            endingpoint=pth[k].b()
            k+=1
            i=startingpoint.position
            while self.border()[i].a()!=endingpoint:
                rr.append(self.border()[i])
                i=(i+1)%len(self.border())
            if k==k0: break
        return GeodesicFigure(rr)
        
    #computes the intersection number between paths a and b given as lists of edges
    def intersection_number(self,a,b):
        BL=len(self.border())  
        rr=0
        lb=[k.label for k in b]
        for k in a:
            if not (k.label in lb or -k.label in lb):
                l=k.costarters()[1:]
                for j in l:
                    if j.label in lb: rr=rr+1/2;  break
                    if -j.label in lb: rr=rr-1/2; break
                l=k.coterminals()[1:]
                for j in l:
                    if j.label in lb: rr=rr+1/2; break
                    if -j.label in lb: rr=rr-1/2; break
        return rr
        
    def represent_loop(self,loop):
        """
        Given a closed path, this function returns a 
        vector representing it in terms of the homology basis
        """
        l=matrix(ZZ,[0 for k in range(len(self.border())//2)])
        for edge in loop:
            ed=edge.label
            l[0,abs(ed)-1]+=ed/abs(ed)
        rr=l*self._helper_matrix
        if rr*self.homology()!=l: print('oops')
        return vector(rr[0])
        
    def g_action(self,g,loop):
        """
        Given a group element and a closed path (as a list of border edges), 
        this function calculates the image of the path under the action of g
        and returns an equivalent path consisting of border edges.
        """
        l=[self.faces[self.face_id[g*e.face.g]]._edges[e.direction] for e in loop]
        return self.equivalent_path(l)
    
        
    def _representation(self,g):
        gloops=[self.g_action(g,loop) for loop in self.loops]
        return matrix([self.represent_loop(loop) for loop in gloops]).transpose()
    
    def __compute_homology(self):
        BL=len(self.border()) 
        #Creates matrix of equations a closed path will have to satisfy
        m=matrix(BL//2,self.number_of_vertices)
        for k in range(BL):
            ed=self.border()[k]
            if ed.label>0:
                m[ed.label-1,ed.a().label]=-1
                m[ed.label-1,ed.b().label]=1
        self.mm=m
        self._homology=m.antitranspose().transpose().kernel().matrix().antitranspose().transpose()
        self._helper_matrix=self._homology.solve_right(identity_matrix(ZZ,self._homology.nrows()))
        self._Genus=self._homology.nrows()//2
        self.loops=[GeodesicFigure(self._homology_basis_element(k)) for k in range(self._homology.nrows())]#basis elements as lists of edges 
        
    def intersection_matrix(self):
        if self._im==0:  #matrix has not been computed
            if self._homology==0: self.__compute_homology() #homology has not been computed
            self._im=matrix(self._homology.nrows())#m_ij will be the intersection number of the ith and jth basis elements
            for k in range(self._homology.nrows()):
                for j in range(k+1,self._homology.nrows()):
                    self._im[k,j]=int(self.intersection_number(self.loops[k],self.loops[j]))
                    self._im[j,k]=-self._im[k,j]
        return self._im  
    
    def __compute_symplectic_group(self):
        self.J,P=self.intersection_matrix().symplectic_form()
        #P is a change of basis which makes the intersection matrix be in standard symplectic form
        Q=P.transpose().inverse()
        symplectic_group_generators=[Q*self._representation(g)*P.transpose() for g in self.generators]
        sgroup=symplectic_group_generators[0].parent()
        #symplectic representations of the given generators
        #here we will create a list with the symplectic representations of all the group elements
        symplectic_group_list=[sgroup(f.word(symplectic_group_generators)) for f in self.faces]
        self.symplectic_group_list=symplectic_group_list
        M=better_basis(symplectic_group_list)#Attempts to find a new basis such that the representations of group elements contain blocks of zeroes.
        self._sgg=[M.transpose().inverse()*k*M.transpose() for k in symplectic_group_generators]
        self._sgl=[M.transpose().inverse()*k*M.transpose() for k in symplectic_group_list]
        self.change_of_basis=M*P

    def symplectic_group_generators(self):
        if self._sgg==0: 
            self.__compute_symplectic_group()
        return self._sgg
    
    def _symplectic_group_list(self):
        if self._sgl==0: 
            self.__compute_symplectic_group()
        return self._sgl

    def relations(self):
        if self._rel==0:
            self._rel=[]
            for k in range(len(self.border())//self.signature[0]):
                g=self.border()[k]
                w=g.face.word*self._ewg[g.direction]*g.rev.face.word.inverse()
                self._rel.append(w)
        return self._rel
                         
    def moebius_invariant(self):
        """
        returns a general matrix solution or a list of equations and a matrix with unknowns
        """
        matrices=self.symplectic_group_generators()
        return moebius_invariant(matrices)

    def moebius_invariant_ideal(self):
        """
        returns a general matrix solution or a list of equations and a matrix with unknowns
        """
        matrices=self.symplectic_group_generators()
        return moebius_invariant_ideal(matrices)

    def _co_starters(self,edgenumber):
        if not self._cs.has_key(edgenumber):
            rr=[]
            next=-self._edges[self._edges.index(edgenumber)-1]
            while next!= edgenumber:
                rr.append(next)
                next=-self._edges[self._edges.index(next)-1]
            self._cs[edgenumber]=rr
        return self._cs[edgenumber]
        
    def _co_terminals(self,edgenumber):
        return [-k for k in self._co_starters(-edgenumber)]
    
    def apply_word_to_face(self,word,f):
        if word.is_one():
            g=self.generators[0]^0
            m=self.fundamental_polygon.rotation_matrices[0]^0
        else:
            g=word(self.generators)
            m=word(self.fundamental_polygon.rotation_matrices)
        f.label=self.face_id[g]
        f.word=word*f.word
        G=g.parent()
        f.g=g*G(f.g)
        f.apply_moebius(m)
        
    def select_faces_by_words(self,L):
        """L should be a list of words in the generators of the group.

        An example of a word is a^3*b^-1*c*a
        
        The result is an object which consists of a 
        face of Poly for each word in the list.
        This object can be shown using the draw() method

        Example:

        P=Poly(SymmetricGroup(4),[4,2,2,2])
        C=P.select_faces_by_words([a,a*b^-1])
        C.draw()
        """
        S=Face(self.fundamental_polygon,self.generators[0]^0,0,self.free_group().0^0)
        picture=GeodesicFigure([deepcopy(S) for k in L])
        for k in range(len(L)):
            self.apply_word_to_face(L[k],picture[k])
        return picture    

    def rational_characters_acting(self):
        if self._rational_characters==None:
            G=self._G
            Ms=self.symplectic_group_generators()[0].parent()
            h1=dict()
            for f in self.faces:
                h1[f.g]=Ms(f.word(self.symplectic_group_generators()))
            pp=self.QG().module_morphism(on_basis=lambda x: h1[x],codomain=Ms)
            cc=G.conjugacy_classes()
            L=[Ms(pp(self.QG()(k.representative()))).trace() for k in cc]
            cf=ClassFunction(G,L)
            cfi=cf.irreducible_constituents()
            T=[vector(k.values()) for k in cfi]
            TQ=[]
            while T!=[]:
                t=T[0]
                GG=t.base_ring().galois_group()
                t_orbit=[t]
                T.remove(t)
                for g in GG:
                    gt=t.apply_map(g)
                    if not gt in t_orbit:
                        t_orbit.append(gt)
                        T.remove(gt)
                c=vector(QQ,sum(t_orbit))
                TQ.append(ClassFunction(G,c))
            self._rational_characters=TQ
        return self._rational_characters
    
    def central_idempotents(self):
        if self._central_idempotents==None:
            G=self.group()
            self._central_idempotents=[]
            T=self.rational_characters_acting()
            Ms=self.symplectic_group_generators()[0].parent()
            h1=dict()
            for f in self.faces:
                h1[f.g]=Ms(f.word(self.symplectic_group_generators()))
            pp=self.QG().module_morphism(on_basis=lambda x: h1[x],codomain=Ms)
            for X in T:
                N=X.irreducible_constituents()[0].degree()
                mX=N/G.order()*sum(QQ(X(g.inverse()))*self.QG()(g) for g in G)
                self._central_idempotents.append(pp(mX))
        return self._central_idempotents

    def primitive_H_idempotent_affording_X(self,X,H=None):
        if H==None:
            H=self._G
        Xi=X.irreducible_constituents()[0]
        e0=Xi.degree()/H.order()*sum(QQ(X(g.inverse()))*self.QG()(g) for g in H)
        if Xi.degree()==1:
            return e0
        L=[k for k in H.conjugacy_classes_subgroups() if k.order()!=H.order()]  #podriamos restringir a subgrupos no-normales?
        Es=[]
        for H in L:
            XH=Xi.restrict(H)
            XHdec=XH.decompose()
            X1=XHdec[0]
            if X1[0]==1:
                rc=rational_character_for_X(X1[1])
                e1=self.primitive_H_idempotent_affording_X(rc,H=H)
                return e1*e0
        raise ValueError("no suitable subgroup found")

def regular_polygon(n,r):
    """returns a geodesic polygon of n sides whose euclidean radius is r"""
    return GeodesicPolygon([r*exp(2*Pi*k/n*i) for k in range(n)])

def rational_character_for_X(X):
    G=X.domain()
    t=vector(X.values()) 
    GG=t.base_ring().galois_group()
    t_orbit=[t]
    for g in GG:
        gt=t.apply_map(g)
        if not gt in t_orbit:
            t_orbit.append(gt)
    c=vector(QQ,sum(t_orbit))
    return ClassFunction(G,c)


def find_generators(G,A):
    return find_generating_vector(G,A)
        
def find_generating_vector(G,A):
    n=len(A)
    As=A[:n-1] #A minus its last entry
    CC=[[k for k in G if k.order()==j]for j in As]
    kk=0
    for k in cartesian_product_iterator(CC): #searches for generators of prescribed orders in the group G
        if prod(k).order()==A[n-1] and G.subgroup([j for j in k]).order()==G.order(): kk=k; break
    if kk==0:
        raise ValueError("Found no generating set")
    l=list(kk)
    l.append(prod(l)^-1)
    return l

def find_all_generators(G,A):
    return find_all_generating_vectors(G,A)

def find_all_generating_vectors(G,A):
    n=len(A)
    As=A[:n-1] #A minus its last entry
    CC=[[k for k in G if k.order()==j]for j in As]
    kk=[]
    for k in cartesian_product_iterator(CC): #searches for generators of prescribed orders in the group G
        if prod(k).order()==A[n-1] and G.subgroup([j for j in k]).order()==G.order(): 
            k1=list(k)
            k1.append(prod(k1)^-1)
            #k1.sort()
            #k2=tuple(k1)
            kk.append(k1)
    if kk==[]:
        raise ValueError("Found no generating set")
    return kk   #  [list(k)for k in kk]

def permutation_representation(G):
    "Finds a permutation representation of the given group acting on cosets"
    g=G._gap_()
    L=g.ConjugacyClassesSubgroups()
    for k in range(len(L)):
        N=g.Core(L[len(L)-k].Representative())
        if N.Order()==1:
            break
    h=L[len(L)-k].Representative()
    N=g.Order()/h.Order()
    G1=PermutationGroup(from_gap_list(SymmetricGroup(N),g.FactorCosetAction(h).Image().GeneratorsSmallest()))
    return G1
    
def perm_rep_gap_group(G):
    return PermutationGroup(gap.IsomorphismPermGroup(G).Image().GeneratorsOfGroup())
    
def third_point(X,Y,d1,d2):
    """
    Given two points X, Y and two hyperbolic distances d1, d2
    it attempts to find a third point Z such that D(X,Z)=d1
    and D(Y,Z)=d2
    """
    if not(abs(X)<1 and abs(Y)<1):
        raise ValueError("The given points are not inside the Poincare disc")
    d3=hyperbolic_distance(X,Y)
    if not(d1+d2>d3 and d1+d3>d2 and d2+d3>d1):
        raise ValueError("Triangle inequality is not satisfied")
    C=GeodesicArc(X,Y)
    C.move(-X)
    alpha=arg(C[1])
    C.rotate(-alpha,CDF(0))
    A2=(RDF(cosh(d1))-1)/2
    A3=(RDF(cosh(d2))-1)/2
    c=real(C[1])
    a=(c*c+A2/(1+A2)-A3/(1+A2)*(1-c*c))/(2*c)
    b=sqrt(A2/(1+A2)-a*a)
    C=GeodesicArc(CDF(0),a+b*i)
    C.rotate(alpha,CDF(0))
    C.move(X)
    return C[1]

def hyperbolic_triangle(d1,d2,d3):
    if not(d1+d2>d3 and d1+d3>d2 and d2+d3>d1):
        raise ValueError("Triangle inequality is not satisfied")
    A1=(RDF(cosh(d1))-1)/2
    c=sqrt(A1/(1+A1))
    X=CDF(0)
    Y=CDF(c)
    Z=third_point(X,Y,d3,d2)
    return GeodesicPolygon([X,Y,Z])

def hyperbolic_distance(X,Y):
    if not(abs(X)<1 and abs(Y)<1):
        raise ValueError("The given points are not inside the Poincare disc")
    return acosh(2*abs(X-Y)^2/((1-abs(X)^2)*(1-abs(Y)^2))+1)    
    
def descparcial(x,n,m,j):
    """
    comienza la descomposicion de x como suma de n fracciones 1/i cuyos denominadores dividen a m y son al menos j (se asume j mayor o igual a 2.
    entrega un par [1/i,x-1/i].
    """
    k=(n/x).floor();
    if n == 1:
        if x == 1/k:
            return [[k]]
        return []
    k0=max((1/x+1).floor(),j)
    d=[i for i in range(k0,k+1) if m.mod(i)==0]
    l=[[x-1/i,i] for i in d]
    return l

def desctotal(x,n,m,j):
    if n == 1:
        return descparcial(x,n,m,j)
    l=descparcial(x,n,m,j)
    r=[]
    for i in l:
        t1=desctotal(i[0],n-1,m,i[1])
        s=[k+[i[1]] for k in t1]
        r.extend(s)
    return r

def desc(x,n,m):
    "Comentario "
    return desctotal(x,n,m,2)
    r=desctotal(x,n,m,2)
    if r!=[]:
        return r	



def riemann_hurwitz(n,g,gamma=None):
    """
    Given a group order n and a genus g gives a list of all compatible signatures.
    This does not verify existence of a group with that signature.
    The optional parameter gamma specifies the genus of the quotient.
    If you have a specific group and you want gamma=0, try suggest_signatures(G,g) instead.
    """
    if gamma==None:
            maxgamma=((g-1)/n+1).floor()
            r=[]
            for gamma in range(maxgamma+1):
                r1=riemann_hurwitz(n,g,gamma=gamma)
                if r1!=[]:
            	    r.append([{'gamma':gamma},r1])    
            return r
    n=ZZ(n)
    k=((g-n*(gamma-1)-1)*4/n).floor()
    k0=(g-n*(gamma-1)-1)*2/n
    if k0 == k0.floor(): k0+=1
    k0=k0.floor()
    rp=[]
    for t in range(k0,k+1):
        c=t+2*(gamma-1)-2*(g-1)/n
        d=desc(c,t,n)
        rp.extend(d)
    return rp


def suggest_signatures(G,g):
    "List of possible signatures for G acting on a surface of genus g with gamma=0"
    r=riemann_hurwitz(G.order(),g,0)
    s=[]
    for i in r:
        if G.is_abelian():
            if G.order()>4*g+4: continue
            if slcm(i)<G.exponent(): continue #a cheap test to discard some signatures on the abelian case
        try: 
            find_generating_vector(G,i)
            s.append(i)
        except(ValueError): 
            pass 
    return s


def slcm(L):
    return min([lcm([L[k] for k in range(len(L)) if k!=j]) for j in range(len(L))])


def XGCD(L):
    if len(L)==1:
        if L[0]>0:
            return L[0],[1]
        if L[0]<0:
            return -L[0],[-1]
    if len(L)==2:
        g,a,b=xgcd(L[0],L[1])
        return g,[a,b]
    g1,A=XGCD(L[1:])
    g,a,b=xgcd(L[0],g1)
    return g,[a]+[b*k for k in A]

def split_space(ListOfMatrices):
    bound=2
    g=ListOfMatrices[0].nrows()//2
    M=matrix([sum([list(k) for k in list(m)],[])for m in ListOfMatrices])
    ListOfMatrices=[matrix(2*g,list(l)) for l in M.row_space().basis()]
    tup=[vector([0 for i in range(2*g)])]
    LW=[]
    P=0
    J=matrix(2*g,2*g)
    for i in range(g):
        J[i,i+g]=1
        J[i+g,i]=-1
    cc=0
    for i in range(1,bound+1):
        tup2=[]
        for l in tup:
            for j in range(2*g):
                if l[j]==0:
                    t=copy(l)
                    t[j]=1
                    tup2.append(t)
        for l in tup2:
            cc+=1
            l0=[k*l.column() for k in ListOfMatrices]
            l1=[(l*J*k)[0] for k in l0]
            if not any(l1): #all entries are 0, so all loops in l0 are disjoint
                m=matrix([k.transpose()[0] for k in l0]).row_space()
                if m.dimension()==g:
                    LW.append(m)
        for j in range(len(LW)):
            for k in range(j+1,len(LW)):
                m=LW[j]+LW[k]
                if m.dimension()==2*g:
                    listo=2
                    P=list(LW[j].basis())
                    P.extend(LW[k].basis())
                    m=matrix(P)
                    if abs(m.det())==1:
                        return m

        if LW!=[]:
            return None
            FullBasis=matrix(2*g,2*g,1).row_space().basis()
            for W in LW:
                m=matrix(2*g)
                for j in range(g):
                    m[j]=W.basis()[j]
                for j in range(g):
                    l=[x for x in FullBasis if m[j]*J*x>0]+[-x for x in FullBasis if m[j]*J*x<0]
                    il=[m[j]*J*x for x in l]
                    il0=list(set(il))
                    il0.sort()
                    l0=[l[il.index(k)] for k in il0]
                    if GCD(il0)!=1: raise ValueError("this shold not happen!!!, tell Antonio")
                    while il0[0]!=1:
                        d,x,y=xgcd(il0[0],il0[1])
                        l0=[x*l0[0]+y*l0[1]]+l0[2:]
                        il0=[d]+il0[2:]
                    m[g+j]=v=l0[0]
                    for k in range(j+1,g+j):
                        m[k]=m[k]-(m[k]*J*v)*m[j]
                if abs(m.det())==1: return m
    return None

def better_basis(ListOfMatrices):
    #return identity_matrix(ListOfMatrices[0].nrows())
    m=split_space(ListOfMatrices)
    g=ListOfMatrices[0].nrows()//2
    J=matrix(2*g,2*g)
    for i in range(g):
        J[i,i+g]=1
        J[i+g,i]=-1
    if m==None: return identity_matrix(ListOfMatrices[0].nrows())
    
    for i in range(g):
        l1=m.submatrix(i,0,g-i,2*g)*J*m.submatrix(g+i,0,g-i,2*g).transpose()
        l2=[GCD(k) for k in l1]
        d,A=XGCD(l2)
        A=matrix(A)
        r=A.nonzero_positions_in_row(0)[0]
        A.swap_columns(0,r)
        m.swap_rows(i,r+i)
        m[i]=sum([A[0][k]*m[i+k] for k in range(A.nrows())])
        l3=m[i]*J*m.submatrix(g+i,0,g-i,2*g).transpose()
        d,A=XGCD(l3)
        if d<0:
            print(l3)
        A=matrix(A)
        r=A.nonzero_positions_in_row(0)[0]
        A.swap_columns(0,r)
        m.swap_rows(g+i,g+r+i)
        m[g+i]=sum([A[0][k]*m[g+i+k] for k in range(A.nrows())])
        for j in range(g):
            if j==i: continue
            a=(m.rows()[g+j]*J*(m.rows()[i].column()))[0]
            b=(m.rows()[j]*J*(m.rows()[g+i].column()))[0]
            m[g+j]=m[g+j]+a*m[g+i]
            m[j]=m[j]-b*m[i]
    return m


def rotation_matrix(alpha,z_0): 
    """
    given a rotation angle alpha and a center z_0 
        (complex number corresponding to the point on the Poincare disc),
    this routine returns a matrix corresponding to the 
    Moebius transformation for the 
    hyperbolic rotation counterclockwise around z_0"""
    m= matrix(2,[exp(i*alpha)-z_0*z_0.conjugate(),z_0*(1-exp(i*alpha)),z_0.conjugate()*(exp(i*alpha)-1),1-z_0*z_0.conjugate()*exp(i*alpha)])
    return m

def moebius_transformation(z,M):
    """Applies the Moebius transformation given by 
    the matrix M to the complex number z"""
    r=(M[0,0]*z+M[0,1])/(M[1,0]*z+M[1,1])
    return r 


def rho(L):
    x=var('__x')
    n=len(L)-1
    def f(x):
        g=sum(asin(sqrt((cos((2*Pi)/L[k])+1)/(sinh(x)^2*(1-cos((2*Pi)/n))+2))) for k in range(n))
        return g+n*asin(sqrt((cos((2*Pi)/n)+1)/(sinh(x)^2*(1-cos((2*Pi)/n))+2)))-Pi/L[n]
    r=find_root(f,0,10)
    return RDF(r)
    
def tesselate_disc(A,bound=0.1,S=None,klein=False,alpha=0,G=SymmetricGroup(1),alphabet="abcdefghijklmnopqrstuvwxyz"):
    if klein: 
        alpha=-Pi/6
    L=FundamentalPolygon(A,beta=alpha,alphabet=alphabet)
    n=len(A)
    letters=L.letters()
    MatSpace=L.rotation_matrices[0].parent()
    #Creates the list of generators adding their inverses
    #There are several representations: group elements, moebius matrices, fuchsian pre-images
    #no entry is created for the last generator, effectively ignoring it
    ewg=sum([[k,k^-1]for k in [letters[k] for k in range(n-1)]],[])
    ng=len(ewg)
    #
    #setup of central part generated by first generator
    generators=find_generating_vector(G,A) if G!=SymmetricGroup(1) else [G.identity()for k in ewg]
    eg=sum([[k,k^-1]for k in generators[:-1]],[])
    g0=generators[0]
    w0=ewg[0]
    N=A[0]
    faces=[Face(L,g0^k,k,w0^k) for k in range(N)]
    face_id=dict([(g0^k,k) for k in range(N)])
    for g in faces:
        g._edges[0].is_border=False
        g._edges[1].is_border=False
    l=N
    ready=False
    while not ready:
        ready=True
        for k in range(ng):
            h=k+(-1 if k%2 else 1)
            for g in range(N):
                if faces[g]._edges[k].is_border==None:
                    w=faces[g].word*ewg[k]
                    cn=faces[g].g*eg[k]
                    cf=Face(L,cn,l,w)
                    try:
                        fl=face_id[cn]
                    except(KeyError):
                        face_id[cn]=l
                        l+=1
                    else:
                        cf.label=fl
                    flag=False
                    for f in faces:
                        if (cf.center()-f.center()).abs()<bound:
                            flag=True
                            f._edges[h].is_border=False
                            faces[g]._edges[k].is_border=False
                            break
                    if flag: continue
                    ready=False
                    faces.append(cf)
                    #faces[g]._edges[k].is_border=False
                    #cf._edges[h].is_border=False
        N=len(faces)
    if S==None:
        picture=GeodesicFigure(faces)
    else:
        picture=[]
        for k in faces:
            c=deepcopy(S)
            m=MatSpace(k.word(L.rotation_matrices))
            c.apply_moebius(m)
            picture.append(c)
        picture=sum(picture[1:],picture[0])
    if klein:
        picture=picture+picture.reflected()    
    return picture
    
def submatrices(M):
    g=M.nrows()//2
    A=M.submatrix(0,0,g,g)
    B=M.submatrix(0,g,g,g)
    C=M.submatrix(g,0,g,g)
    D=M.submatrix(g,g,g,g)
    return A,B,C,D
    
def moebius_invariant(L,Blocks=1,polarization=[]):
    """
    Given a list L of symplectic matrices this function returns a polynomial ring R, and ideal I and a matrix Z with coefficients in R whose image in R/I is invariant under Moebius action by the given matrices
    """
    g=L[0].nrows()//2
    if polarization==[]:
        polarization=[1 for _ in range(g)]
    polarization+=[1 for _ in range(g)]
    alfa=diagonal_matrix(polarization)
    matrices=[alfa*k*alfa.inverse() for k in L]
    n=len(matrices)
    a=[matrix([k[:g]for k in m[:g]])for m in matrices]
    b=[matrix([k[g:2*g]for k in m[:g]])for m in matrices]
    c=[matrix([k[:g]for k in m[g:2*g]])for m in matrices]
    d=[matrix([k[g:2*g]for k in m[g:2*g]])for m in matrices]
    gd=g//Blocks
    Nvars=Blocks*gd*(gd+1)//2
    if Nvars>15:
        raise ValueError("too many variables to be realistic")
    R=PolynomialRing(QQ,'x',Nvars,order='degrevlex')
    l=[R.gen(k) for k in range(Nvars)]
    l.reverse()
    gens=copy(l)
    Zs=[]
    for ll in range(Blocks):
        Z1=matrix(R,gd)
        for k in range(gd):
            for j in range(k,gd):
                Z1[k,j]=Z1[j,k]=l.pop()
        Zs.append(Z1)
    Z=block_diagonal_matrix(Zs)               
    l=[(Z*d[k]+b[k])-(Z*c[k]+a[k])*Z for k in range(n)]
    ecs=list(set([mm for _ in l for mm in _.list() ]))
    if forall([k.degree() for k in ecs],lambda x: x<=1)[0]:
        ecs=[SR(k) for k in ecs]
        l=[SR(k) for k in gens]
        if SR(0) in ecs:
            ecs.remove(SR(0))
        s=solve(ecs,l,solution_dict=true)
        Zs=[]
        for ll in range(Blocks):
            Z1=matrix(SR,gd)
            for k in range(gd):
                for j in range(k,gd):
                    Z1[k,j]=Z1[j,k]=l.pop()
            Zs.append(Z1)
        Z=block_diagonal_matrix(Zs)
        return [Z.substitute(s[0])]
    else:
        I=R.ideal(ecs)
        l=[SR(k) for k in gens]
        Zs=[]
        for ll in range(Blocks):
            Z1=matrix(R,gd)
            for k in range(gd):
                for j in range(k,gd):
                    Z1[k,j]=Z1[j,k]=l.pop()
            Zs.append(Z1)
        Z=block_diagonal_matrix(Zs)
        #print "Since equations are not linear, it is not known how hard it can be to solve them.  We return an ideal containing all the necessary equations (as polynomials) and a matrix of with the form of the invariant riemann matrices"
        return I,Z

def moebius_invariant_ideal(L):
    """
    Given a list L of symplectic matrices this function returns a polynomial ring R, and ideal I and a matrix Z with coefficients in R whose image in R/I is invariant under Moebius action by the given matrices
    """
    matrices=L
    g=L[0].nrows()//2
    n=len(matrices)
    a=[matrix([k[:g]for k in m[:g]])for m in matrices]
    b=[matrix([k[g:2*g]for k in m[:g]])for m in matrices]
    c=[matrix([k[:g]for k in m[g:2*g]])for m in matrices]
    d=[matrix([k[g:2*g]for k in m[g:2*g]])for m in matrices]
    R=PolynomialRing(QQ,'x',(g*(g+1))// 2,order='lex')
    l=[R.gen(k) for k in range(g*(g+1)//2)]
    l.reverse()
    gens=copy(l)
    Z=matrix(R,g)
    for k in range(g):
        for j in range(k,g):
            Z[k,j]=Z[j,k]=l.pop()               
    l=[(Z*d[k]+b[k])-(Z*c[k]+a[k])*Z for k in range(n)]
    ecs=list(set([mm for _ in l for mm in _.list() ]))
    I=R.ideal(ecs)
    return I
    

def SmallGroup(a,b):
    """
    Wrapper for GAP function which returns a permutation group from the small groups library (the small groups library has to be installed in GAP)
    gap.eval('g:=Image(IsomorphismPermGroup(SmallGroup(%d,%d)))'%(a,b))
    return PermutationGroup(gap_group='g')
    """
    g=libgap.SmallGroup(a,b)
    return g

def QuotientHomomorphism(G,N):
    ghom=G._gap_().FactorCosetAction(N._gap_())
    quot=PermutationGroup(gap_group=ghom.Image())
    h=PermutationGroupMorphism_im_gens(G, quot,map(ghom.Image,G.gens()))
    return h

from sage.groups.perm_gps.permgroup import from_gap_list
def AutOrbit(L,G,inner=False):
    """
    Given a list of elements generating a group G, the function lets the automorphism group of G act on the list producing a new list of elements of G.  Finally it returns the orbit of L under this action.
    """
    gap.eval('G := Group(' + str(L) + ')')
    gap.eval('aut := AutomorphismGroup(G)')
    if inner:
        gap.eval('aut := InnerAutomorphismsAutomorphismGroup(aut)')
    gap.eval('L:=List(aut,y->List('+str(L)+',x->x^y))')
    max=eval(gap.eval('Length(L)'))
    LL=[gap.eval('Print(L['+str(k)+'])') for k in range(1,max+1)]
    #G=PermutationGroup(L)
    LL3=[from_gap_list(G,k.replace('\n','').replace(' ','')) for k in LL]
    M=[[G(k) for k in m]for m in LL3]
    #for k in M:
    #    k.sort()
    M1=[tuple(k) for k in M]
    M2=list(set(M1))
    M2.sort()
    return M2

def Representatives(L,G,inner=False):
    """
    given a list L=[v_1,v_2,...v_k] whose elements are lists or tuples of group elements, it computes the orbits of the action of Aut(G)XB 
    on the set of v's and returns a representative of each orbit.  It doesn't matter if the orbit of some of the v's is not completely contained 
    in L to begin with.
    """ 
    M=set([])
    C=[]
    for k in range(len(L)):
        if not tuple(L[k]) in M:
            C.append(L[k])
            M1=set(Orbit(L[k],G,inner=inner))
            M=M.union(M1)
    return C

def q(k,v):
    """
    Performs the kth braid operation on the vector v of group elements.  Does not change v.
    """
    l=list(v[:k])
    c=[v[k]*v[k+1]*v[k]^-1,v[k]]
    r=list(v[k+2:])
    return l+c+r

def BraidOrbit(X):  #X is a list of lists of group elements
    terminamos=False
    L=set([tuple(k) for k in X])
    r=len(X[0])-1
    while terminamos==False:
        m=len(L)
        M=copy(L)
        for k in range(r):
            for w in L:
                M.add(tuple(q(k,w)))
        L=M
        if len(L)==m:
            terminamos=True
    return list(L)

def Orbit(L,G,inner=False):
    return BraidOrbit(AutOrbit(L,G,inner=inner))

def find_generator_representatives(G,A,inner=False):
    n=len(A)
    As=A[:n-1] #A minus its last entry
    CC=[[k for k in G if k.order()==j]for j in As]
    kk=[]
    for k in cartesian_product_iterator(CC): #searches for generators of prescribed orders in the group G
        if prod(k).order()==A[n-1]: #checks that the product is one
            k1=list(k)
            k1.append(prod(k1)^-1)
            kk.append(k1)
    if kk==[]:
        raise ValueError("Found no generating set")
    kk1=Representatives(kk,G,inner=inner)  #reduces the list by the action of automorphisms of G and the Braid group.
    kk2=[]
    for k in kk1:
        if G.subgroup([j for j in k]).order()==G.order():  #removes those tuples which do not generate the group G.
            kk2.append(k)
    if kk2==[]:
        raise ValueError("Found no generating set")        
    return kk2  


def positive_Tietze(T,orders):
    S=[]
    X=list(T)
    while X!=[]:
        a=X.pop(0)
        k=1
        while X!=[] and X[0]==a:
            k+=1
            X.pop(0)
        if a>0:
            S.extend([a]*k)
        else:
            S.extend([-a]*(orders[-a-1]-k))
    return tuple(S)

def as_word(x,G,gens=None,positive=True):
    if gens==None:
       gens=G.gens()
    orders=[k.order() for k in gens]
    Gfp=G.as_finitely_presented_group()
    iso = libgap.IsomorphismFpGroupByGenerators(G, gens)
    u=Gfp(libgap.Image(iso,x)).Tietze()
    if positive:
    	w=positive_Tietze(u,orders)
    return Gfp(w)

def find_generator_representatives_as_words(G,A,inner=False):
    gens=find_generator_representatives(G,A,inner=inner)
    return [[as_word(k,G) for k in l] for l in gens]



def char_poly_from_traces(pp,x):
    dd=len(pp)-1
    s=[1]
    for k in range(1,dd+1):
        r=-1/k*sum((-1)^l * s[k-l]*pp[l] for l in range(1,k+1))
        s.append(r)
    ans=sum(s[k]*(-1)^k*x^(dd-k)for k in range(dd+1))
    return ans

      

def n11(g,chi,z,cc):
    G=g.parent()
    K=z.parent()
    t=z^(z.multiplicative_order()//g.order())
    p=[K(chi(g^k)) for k in range(chi.degree()+1)]
    S.<x>=K[]
    px=char_poly_from_traces(p,x)
    rans=[]
    for alpha in range(g.order()):
        f=copy(px)
        ans=0
        while f(x=t^alpha)==0:
            f=f.derivative()
            ans+=1
        rans.append(ans)
    return rans
    
def CW(G,v):
    K.<z>=CyclotomicField(G.exponent())
    T=G.irreducible_characters()
    cc=G.conjugacy_classes()
    e=[k.order() for k in v]
    K=z.parent()
    ans=0*T[0]
    for i in range(1,len(T)):
        mu=-T[i].degree()
        for j in range(len(v)):
            nn11=n11(v[j],T[i],z,cc)
            for alpha in range(1,e[j]):
                mu+=((1-alpha/e[j])*nn11[alpha])
        ans+=mu*T[i]
    return ans


def SerreFormula(G,X):
    L=dict([(g,(libgap(g)^X).sage())for g in G])
    N=1/(2*G.order())*sum([L[g]^2+L[g^2] for g in G])
    return N

def Jmatrix(g=1,D=None):
    if D==None:
        D=[1]*g
    g=len(D)
    m1=diagonal_matrix(D)
    m0=zero_matrix(g)
    m=block_matrix(2,[m0,m1,-m1,m0],subdivide=False)
    return m

def moebius_action(DZ,M):
    DZM=DZ*M
    d=DZ.nrows()
    A=DZ.submatrix(0,0,d,d)
    B=DZM.submatrix(0,0,d,d)
    return A*B.inverse()*DZM

def reduce_signature(L):
    S=list(set(L))
    S.sort()
    S.reverse()
    T=[]
    for k in S:
        if L.count(k)==1:
            T.append(k)
        elif L.count(k)==2:
            T.append(k)
            T.append(k)
        else:
            T.append(str(k)+"^"+str(L.count(k)))
    return T


def intermediate_covers(v,G,H,gwG,verbose="English"):
    if v[0] in G:
        l=[G.subgroup([k])for k in v]
    else:
        l=v
    a=len(l)
    ng=[G.normalizer(k) for k in l]
    rtng=[[k[0] for k in G.cosets(g)]for g in ng]
    s=[[l[k].conjugate(g) for g in rtng[k]]for k in range(a)]
    t=[[H.intersection(k)for k in j]for j in s]
    index=[ng[k].order()/l[k].order() for k in range(a)]
    B=[sum(index[k]*(j.order()-1)for j in t[k]) for k in range(a)]
    b=sum(B)
    bi=sum(1-1/k.order() for k in l)
    gw=G.order()*(gwG-1)+1+(G.order()*bi)/2
    gwH=(gw-1-b/2)/H.order() +1
    u=[[index[k]*j.order() for j in t[k]]for k in range(a)]
    u2=[[[j,k.count(j)]for j in set(k)]for k in u]
    sig=[[[j[0]/index[k],j[0]*j[1]/H.order()]for j in u2[k]]for k in range(a)]
    sig2=[[[j[0]/index[k]]*int(j[0]*j[1]/H.order())for j in u2[k]]for k in range(a)]
    m=[[]for k in range(a)]
    for k in range(a):
        for j in sig[k]:
            m[k].extend([l[k].order()/j[0]for _ in range(j[1])])
    core=PermutationGroup(gap_group=G._gap_().Core(H))
    if verbose=="Spanish":
        print
        print("genero de W =",gw)
        print("La signatura de W/H es:")
        print("genero de W/H: ",gwH)
        for i in range(a):
            print(sig[i]," que estan sobre el punto marcado por l[%d]"%i)
        print
        print("Para el cubrimiento W/H-->W/G la estructura de ciclos es:")
        for i in range(a):
            print("l[%d]"%i,"-->",m[i])
        print("El indice de H en G es: ",G.order()/H.order())
        print("Core(H)= ",core.id())
    if verbose=="English":
        print
        print("genus of W =",gw)
        print("The signature of W/H is:")
        print("genus of W/H:",gwH)
        for i in range(a):
            print(sig[i]," that are on the point marked by l[%d]"%i)
        print
        print("For the cover W/H-->W/G the structure of the cycles is:")
        for i in range(a):
            print("l[%d]"%i,"-->",m[i])
        print("The index of H in G is: ",G.order()/H.order())
        print("Core(H)= ",core.id())
    s2=flatten(sig2)
    s2=[k for k in s2 if k>1]
    s2.sort()
    s2.reverse()
    return s2

def double_transversal(G,H):
    L=set(G)
    S=[[G.identity()]]
    L=L.difference(set(H))
    while len(L)>0:
        x=L.pop()
        D=set([G(k.sage()) for k in libgap.DoubleCoset(H,x,H).List()])
        L=L.difference(D)
        S1=[]
        while len(D)>0:
            y=D.pop()
            S1.append(y)
            for h in H:
                D.discard(y*h)
                D.discard(h*y)
        S.append(S1)
    return S
    
def kbword_reducer(F):
    G=F.free_group()
    ng=G.ngens()
    gens=G.gens()
    A=F.signature
    rels=[gens[k]^A[k] for k in range(ng)]+[prod(gens)]
    Gq=G.quotient(rels)
    K=Gq.gap().KBMAGRewritingSystem()
    K.AutomaticStructure()
    def f(w):
        return G(K.ReducedWord(w.gap()))
    return f

def braid_max_to_front(v):
    """
    Given a vector v of group elements, braid operations are applied to bring an element of maximal order to the front.
    """
    v0=[k.order() for k in v]
    m=max(v0)
    l=v0.index(m)
    w=copy(v)
    for k in range(l):
        w=q(l-k-1,w)
    return w

def braid_sort(v):
    """
    Given a vector v of group elements, braid operations are applied to sort it so orders are decreasing.
    """
    if len(v)<=1:
        return v
    w=braid_max_to_front(v)
    return [w[0]]+braid_sort(w[1:])


def canonical_generating_vector(v):
    """
    Given a group H and a generating vector, returns a canonical generating vector in a "SmallGroup" isomorphic to H
    """
    H=PermutationGroup(v)
    Hid=H.id()
    H0=SmallGroup(Hid[0],Hid[1])
    f=H.isomorphism_to(H0)
    L=find_generator_representatives(H0,[k.order() for k in v])
    w=tuple([f(k) for k in v])
    for v0 in L:
        if w in Orbit(v0,H0):
            return v0
            
def canonical_generating_vector_as_words(v):
    gens=canonical_generating_vector(v)
    H0=gens[0].parent()
    return [as_word(k,H0) for k in  gens]

def AccionAutomorfismos(L,G):
    M=set([])
    C=[]
    for k in range(len(L)):
        if not tuple(L[k]) in M:
            C.append(L[k])
            M1=set(AutOrbit(L[k],G))
            M=M.union(M1)
    return C
