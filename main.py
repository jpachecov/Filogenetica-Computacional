 #!/usr/bin/python
 # -*- coding: utf-8 -*-
 # author: Jean Pierre Pacheco Avila (jpachecov)

from Bio import SeqIO
import numpy as np
from Bio import Phylo
from cStringIO import StringIO

homo = ''
nomascus = ''
paniscus = ''
trogloditas = ''
pongo = ''

for homo in SeqIO.parse("homo_sapiens", "fasta"):
    print(homo.id)
    #print(repr(homo.seq))
    #print(len(homo))

for nomascus in SeqIO.parse("nomascus_leucogenys", "fasta"):
    print(nomascus.id)
    #print(repr(nomascus.seq))
    #print(len(nomascus))

for paniscus in SeqIO.parse("pan_paniscus", "fasta"):
    print(paniscus.id)
    #print(repr(paniscus.seq))
    #print(len(paniscus))

for trogloditas in SeqIO.parse("pan_troglodytes", "fasta"):
    print(trogloditas.id)
    #print(repr(trogloditas.seq))
    #print(len(trogloditas))

for pongo in SeqIO.parse("pongo_abelii", "fasta"):
    print(pongo.id)
    #print(repr(pongo.seq))
    #print(len(pongo))

def isPurine(base):
	return base in ['a','g']

def isPyrimidine(base):
	return base in ['c','t']


def isTransition(base_1, base_2):
	if (isPurine(base_1) and isPurine(base_2)) or (isPyrimidine(base_1) and isPyrimidine(base_2)):
		return True
	else :
		return False 

def isTransversion(base_1, base_2):
	if (isPurine(base_1) and isPyrimidine(base_2)) or (isPyrimidine(base_1) and isPurine(base_2)):
		return True
	else :
		return False 


# Obtiene la distancia genetica de dos secuencias usando K80
# Suponemos que ambas secuencias tienen la misma longitud
def get_genetic_distance(seq1, seq2):

	if len(seq1) != len(seq2): return

	num_transition = 0
	num_transversion = 0
	cuantos = len(seq1)

	for i in range(0, len(seq1)):
		base_1 = seq1[i]
		base_2 = seq2[i]
	
		if base_1 != base_2:
	
			if isTransition(base_1, base_2): 
				num_transition += 1.
			else :
				if isTransversion(base_1, base_2):
					num_transversion += 1.

	p = num_transition / cuantos
	q = num_transversion / cuantos

	return -0.5 * np.log((1 - 2*p - q) * np.sqrt(1 - 2*q))


# Obyiene la matriz de distancias geneticas a partir de una
# lista de secuencias genomicas
def get_genetic_distance_matrix(secuencias):
	M =  np.zeros((len(secuencias), len(secuencias)))

	for i in range(0, len(secuencias)):
		for j in range(0, len(secuencias)):
			if i < j:
				M[i,j] = get_genetic_distance(secuencias[i], secuencias[j])
	return M


# Una clase que representa al algoritmo Neighbor-Joininig
# Se hace uso de una estructura arbol para representar a los arboles filogeneticos
# qu originalmnte tienen forma de estrella y al final del algoritmo todos los subarboles
# que cuelgan de la raiz son binarios.

class NeighborJoining:

	def __init__(self, matriz, arbol, secuencias):
		self.D = matriz
		self.Q = matriz
		self.T = arbol


		# Para el nombre
		self.S = secuencias

		# Nodos que el algoritmo esta considerando actualmente
		self.currentNodes = arbol.getHijos()
		self.crtNodeNames = map(lambda x: x.nombre, self.currentNodes)
		#self.llenaArbol()

	# Dado un arbol filogenetico inicial en forma de estrella
	#def llenaArbol():

	# Implementacion del algoritmo
	def start(self):
		

		# Las hoas restastantes al principio son todas
		hojas_restantes = self.T.hijos
		
		D = self.D
		# Mientras podamos seguir uniendo
		self.iteracion = 1
		while(len(hojas_restantes) > 2):

			self.update_Q()

			#print "hojas res", map(lambda x: x[1].nombre, hojas_restantes)
			
			# Obtenemos los indices del valor minimo
			minimo = self.getMinimo()

			# Obtenemos los nodos a unir
			node_1 = self.T.getNode(hojas_restantes[minimo[0]].nombre)
			node_2 = self.T.getNode(hojas_restantes[minimo[1]].nombre)

			# Unimos nodos dentro de arbol, esto hara que las aristas
			# al origen disminuyan en 1

			nuevo = self.T.joinNodes(node_1, node_2)

			# obtnemos los hijos del nodo nuevo con sus distancias
			# (d, Nodo)
			d_node1 = nuevo.hijos[0]
			d_node2 = nuevo.hijos[1]
			
			# Actualizamos sus distancias hacia el nuevo nodo
			# Se puede hcar mas exacto

			d_node1.valor = 0.5 * D[minimo[0],minimo[1]]
			d_node2.valor = 0.5 * D[minimo[0],minimo[1]]

			# Se vuelve a calcular una nueva matriz de distancias
			# tniendo en cuenta el nuevo nodo
			# Actuqalizamos matriz D
			
			self.update_D(nuevo, minimo)

			#Actualizamos el conjunto de hojas que aun podemos unir
			hojas_restantes = self.T.hijos

			self.iteracion += 1

		#print "hojas res", map(lambda x: x[1].nombre, hojas_restantes)


	def update_Q(self):
		print 'Q', self.Q
		self.Q = self.D
		"""""
		def suma(i,k):
			sum_ = 0
			for s in range(i, k):
				sum_ += self.D[i,s]
			return sum_

		n = len(self.Q)
		for i in range(0, n):
			for j in range(0, n):
				if i < j:
					self.Q[i,j] = (n - 2) * self.D[i,j] - suma(i, n) - suma(j, n)
		"""""

	# Recibe el nuevo nodo y los incices de los nodos que se unieron
	# Actualiza la matriz D de distancias

	def update_D(self, nuevo, indices):

		# Si A es una matriz de numpy A[i] es un renglon de la matriz A

		"""""
		print "nuevo", nuevo
		viejos_nodos = self.currentNodes
		viejos_nombres = self.crtNodeNames
		distancias = []

		print "ALGO ", nuevo
		print viejos_nombres

		"""""

		print len(self.T.hijos)
		print self.T.hijos
		"""""
		for i in range(0, len(self.S)):
			for j in range(0, len(self.S)):
				if i < j and i != i_ and j != j_ :
					pass
					#distancias.append({'i':})
		"""""

		# Actualizamos valores con el nvo nodo
		n = len(self.D) - 1
		D = np.zeros((n, n))
		
		#print D

		return 0


	# En la primera iteracion obtiene el indice del valor minimo en la matriz Q
	# En las demas iteraciones esta busqueda se hace revisando los nodos 
	# que cuelgan de la raiz del arbol filognetico.

	# Devuelve el par de nodos 

	def getMinimo(self):
		minim = 1000000000
		min_i = 0
		min_j = 0
		for i in range(0, self.Q.shape[0]):
			for j in range(0, self.Q.shape[0]):
				if i < j:
					if Q[i,j] < minim:
						min_i = i
						min_j = j
						minim = Q[i,j]
		
		return (min_i, min_j)

"""""
		if self.iteracion == 1:

			minim = 1000000000
			min_i = 0
			min_j = 0
			for i in range(0, self.Q.shape[0]):
				for j in range(0, self.Q.shape[0]):
					if i < j:
						if Q[i,j] < minim:
							min_i = i
							min_j = j
							minim = Q[i,j]
			
			return (min_i, min_j)
		else:
			anterior = ''
			minimo = self.T.hijos[0]
			for i in range(0, len(self.T.hijos)):
				if self.T.hijos[i].valor < minimo.valor:
					anterior = minimo
					minimo = self.T.hijos[i]
			return (minimo, anterior)
"""""
class Arista:
	def __init__(self):
		self.peso = 0
	def getPeso(self):
		return self.peso
	def setPeso(nuevo):
		self.peso = nuevo

# En este caso el valor sera usado como distancia hacia el nodo padre
class Nodo:
	def __init__(self, nombre, valor):
		self.nombre = nombre
		self.valor = valor
	def getNombre(self):
		return self.nombre
	def getValor(self):
		return self.valor
	def tostring_1(self):
		print self.nombre
		return '' + str(self.valor) + ':(' + self.nombre + ':)'

	def tostring_2(self):
		print self.nombre
		return '' + self.nombre + ':' + str(self.valor) + ':'



# Un arbol es vacio o un nodo con muchos arboles como hijos
# Este arbol inicialmente es una estrella y el algoriitmo NJ lo
# va modificando durante su ejecucion.
# AL final de la ejecucuion podemos llamar a tostring() para
# obtener un recorrido DFS sobre el arbol obteniendo un formato newick



class Arbol_NJ:
	def __init__(self, esRaiz, valor):
		# Una lista de objetos de tipo tuple (d, nodo)
		self.hijos = []
		if(esRaiz):
			self.esRaiz = True
			self.raiz = Nodo('Origen', valor)

		else:
			self.esRaiz = False	
			self.nodo = self.raiz = Nodo('O_', valor)
		#self.aristas
		self.node_index = 1
		self.valor = valor
	def init(self, nombres):
		for name in nombres:
			self.hijos.append(Nodo(name, -1))
		return 0

	def getNombre(self):
		return self.raiz.getNombre()

	def getDistanciaAlPadre(self):
		if self.esRaiz :
			return -1
		else:
			return self.nodo.valor

	# Aplica un racorido en prufundidad (DFS) en el arbol
	# para obtener una forma aceptada por byopython
	# @@@@ ESTE ES UN EJEMPLO DE UNA CADENA BIEN CONSTRUIDA PARA QUE LA IMAGEN
	# SE MUSTRE COMO EN LA DE LA PRACTICA 
	# (2:(A:),5:(2:(5:(B:),3:(C:)),4:(D:)))


	#(A:0.1,B:0.2,(C:0.3,D:0.4):0.5); 
	def tostring_1(self):
		if self.esHoja():
			print "es hoja"
			# es una hoja
			return self.raiz.tostring_1()

		else:
			print "nodo"
			# es un nodo
			cadena = ''
			if self.esRaiz:
				pass
			else:
				cadena += str(self.nodo.valor) + ":"
			cadena
			cadena += '('
			for hijo in self.hijos:
				cadena += hijo.tostring_1() + ','
			cadena = cadena[:-1]
			cadena += ')'
			return cadena
	


	def tostring_2(self):
		if self.esHoja():
			print "es hoja"
			# es una hoja
			return self.raiz.tostring_2()

		else:
			print "nodo"
			# es un nodo
			cadena = ''
			if self.esRaiz:
				cadena += self.raiz.nombre + ':' + str(self.raiz.valor)+":"
				pass
			else:
				cadena += (self.nodo.nombre + ":" + str(self.nodo.valor) + ":")
			
			cadena += '('
			for hijo in self.hijos:
				cadena += hijo.tostring_2() + ','
			cadena = cadena[:-1]
			cadena += ')'
			#
			#O:0.1:(A:0.2:(B:0.3:,C:0.4:), D:0.1:)
			return cadena
	
	def getHijos(self):
		return self.hijos

	def esHoja(self):
		return len(self.hijos) == 0
	def addNode(self, node):
		self.hijos.append(node)
		#print "----------------"
		#print node

		#node[1].nombre += str(self.node_index)
		self.node_index += 1
	def removeNode(self,node):
		self.hijos.remove(node)

	
	def getNode(self, id):
		for node in self.hijos:
			if node.getNombre() == id:
				return node
		return None

	def removeNode(self, node):
		self.hijos.remove(node)
	# @param actual nuevo
	def replace(self, actual, nuevo):
		self.hijos[self.hijos.index(actual)] = nuevo

	# REVISAAAR ESTE!!!!!!!
	def joinNodes(self, nodo1, nodo2):
		
		#creamos nuevo nodo con hojas
		nuevo = Arbol_NJ(False, -1)
		nuevo.addNode(nodo1)
		nuevo.addNode(nodo2)
		self.replace(nodo1, nuevo)
		self.removeNode(nodo2)


		return nuevo

	def aplicaRecorrido(self, cual):
		if cual == 'newick_1':
			return self.tostring_1()
		if cual == 'newick_2':
			return self.tostring_2()




"""""
# DESCOMENTAR PARA PROBAR LA ESTRUCTURA Y MOSTRAR IMAGENES DE ABBOLES filogeneticos
# 

# Creamos un arbol 

padre = Arbol_NJ(True, 0.1)
padre.addNode(Nodo('A', 2.))
otro = Arbol_NJ(False, 5.)
otro_2 = Arbol_NJ(False, 2.)
otro_2.addNode(Nodo('B', 5.))
otro_2.addNode(Nodo('C', 3.))

otro.addNode(otro_2)

otro.addNode(Nodo('D', 4))

padre.addNode(otro)

treedata = padre.aplicaRecorrido('newick_1')
treedata_2 = padre.aplicaRecorrido('newick_2')

handle = StringIO(treedata)
tree = Phylo.read(handle, "newick")
Phylo.draw(tree)

handle_2 = StringIO(treedata_2)
tree_2 = Phylo.read(handle_2, "newick")
Phylo.draw(tree_2)
Phylo.show()
"""""

# Asi se iniciaria


# Secuncias leidas de archivos
secs = [homo, nomascus, paniscus, trogloditas, pongo]

# Etiquetas para identificar las secuencias iniiales
ets = ['A','B','C','D','E']
# Correspondencias
corresp = {}
# Lenado de correspondecisa
for i in range(0,len(ets)):
	corresp[i] = secs[i]


# Obtenemos la matriz inicial de distancias geneticas
# a partir de las secuecnias leidas de los archivos
Q = get_genetic_distance_matrix(secs)

# Raiz del Arbol filogenetico inicial
padre = Arbol_NJ(True, 0.1)
# lo llenamos
for key in corresp:
	padre.addNode(Nodo(key, -1))

# iniciamos instancia de algoritmo
algoritmo_NJ = NeighborJoining(Q, padre, secs)

# Lo eejecutamos
algoritmo_NJ.start()



