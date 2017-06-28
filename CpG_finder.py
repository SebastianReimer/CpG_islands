#!/usr/bin/python



import re #modul fuer regulaere Ausdruecke
from math import log10 #10er Logarithmus-Modul
import matplotlib.pyplot as plt #Visualisierung



###############
###Functions###
###############

###function um Sequenzen von fasta-Datei in Array zu speichern###
#function um eine Sequenz Buchstabe fuer Buchstabe aus einem FASTA-Format in einem Array seq zu speichern.
#moegliche Probleme: Falls Leerzeilen in fasta-file --->Programm gibt Fehler aus
def fasta2array(fastafile):
	i = -1
	seq = []
	with open(fastafile) as f:
		for line in f:
	#		line = re.sub("\s", "", line)	#Entfernt tabs, whitespace nich unbedingt notwendig
			line = line.strip()				#entfernt Zeilenumbruch
			r = list(line)					#einzelne Elemente der Zeile werden Zeichen fuer Zeichen aufgetrennt
			if r[0] == ">":					#prueft ob erstes Zeichen ">"
				i += 1
				seq.append("")				#neuer Leerstring wird hinzugefuegt
			else:
				seq[i] += line				#komplette Proteinsequenz wird konkateniert
		f.close()							#Filepointer schliessen	
	return seq								#nur erste eingelesen Sequemz wird zueckgegeben

###function, um Anzahl der Nukleotide mit deren Nachfolgenukleotiden zu zaehlen; nt1: erste Nukleotid, nt2: Folgenukleotid###
	
def count_nt(array, nt1, nt2):	
	i = 0
	count = 0
	while i < len(array)-1:
		if array[i] == nt1 and array[i+1] == nt2:
			count += 1
		i += 1		
	return count


##############################################
###Aufgabe 1: Transitionsmatrizen a+ und a-###
##############################################

### a(+)-Modell###
##################

cpg = list(fasta2array("CpG.fasta")[0])	#erste Array-element buchstabenweise gesplittet

#Anzahl berechnen von Nuktleotid1 zu Nachfolgenukleotid2 und Speicherung im Transitionsdictionary a_plus_anzahl
a_plus_anzahl = {}
a_plus_anzahl["A"] = {	"A": count_nt(cpg, "A", "A"),
						"C": count_nt(cpg, "A", "C"),
						"G": count_nt(cpg, "A", "G"),
						"T": count_nt(cpg, "A", "T"),	}

a_plus_anzahl["C"] = {	"A": count_nt(cpg, "C", "A"),
						"C": count_nt(cpg, "C", "C"),
						"G": count_nt(cpg, "C", "G"),
						"T": count_nt(cpg, "C", "T"),	}

a_plus_anzahl["G"] = {	"A": count_nt(cpg, "G", "A"),
						"C": count_nt(cpg, "G", "C"),
						"G": count_nt(cpg, "G", "G"),
						"T": count_nt(cpg, "G", "T"),	}

a_plus_anzahl["T"] = {	"A": count_nt(cpg, "T", "A"),
						"C": count_nt(cpg, "T", "C"),
						"G": count_nt(cpg, "T", "G"),
						"T": count_nt(cpg, "T", "T"),	}

# Summe aus Anzahl aller vorkommenden Nukleotide c+
sum_count_c_plus = 0;

for x in a_plus_anzahl:
	for y in a_plus_anzahl[x]:
		sum_count_c_plus += a_plus_anzahl[x][y]

sum_count_c_plus = float(sum_count_c_plus)	#Typkonversion, um bei Berechungen genauere Zahlen zu erhalten


# Transitionsmatrix in Form von Dictionary a_plus berechnen
a_plus = {}

a_plus["A"] = {	"A": count_nt(cpg, "A", "A") / sum_count_c_plus,
				"C": count_nt(cpg, "A", "C") / sum_count_c_plus,
				"G": count_nt(cpg, "A", "G") / sum_count_c_plus,
				"T": count_nt(cpg, "A", "T") / sum_count_c_plus, }

a_plus["C"] = {	"A": count_nt(cpg, "C", "A") / sum_count_c_plus,
				"C": count_nt(cpg, "C", "C") / sum_count_c_plus,
				"G": count_nt(cpg, "C", "G") / sum_count_c_plus,
				"T": count_nt(cpg, "C", "T") / sum_count_c_plus, }

a_plus["G"] = {	"A": count_nt(cpg, "G", "A") / sum_count_c_plus,
				"C": count_nt(cpg, "G", "C") / sum_count_c_plus,
				"G": count_nt(cpg, "G", "G") / sum_count_c_plus,
				"T": count_nt(cpg, "G", "T") / sum_count_c_plus, }

a_plus["T"] = {	"A": count_nt(cpg, "T", "A") / sum_count_c_plus,
				"C": count_nt(cpg, "T", "C") / sum_count_c_plus,
				"G": count_nt(cpg, "T", "G") / sum_count_c_plus,
				"T": count_nt(cpg, "T", "T") / sum_count_c_plus, }

print("Transitionsmatrix a(+)-Modell")
print(a_plus)
print("\n")

### a(-)-Modell###
##################

bckgrnd = list(fasta2array("background.fasta")[0])	#erste Array-element buchstabenweise gesplittet

#Anzahl berechnen von Nuktleotid1 zu Nachfolgenukleotid2 und Speicherung im Transitionsdictionary a_minus_anzahl
a_minus_anzahl = {}
a_minus_anzahl["A"] = {	"A": count_nt(bckgrnd, "A", "A"),
						"C": count_nt(bckgrnd, "A", "C"),
						"G": count_nt(bckgrnd, "A", "G"),
						"T": count_nt(bckgrnd, "A", "T"),	}

a_minus_anzahl["C"] = {	"A": count_nt(bckgrnd, "C", "A"),
						"C": count_nt(bckgrnd, "C", "C"),
						"G": count_nt(bckgrnd, "C", "G"),
						"T": count_nt(bckgrnd, "C", "T"),	}

a_minus_anzahl["G"] = {	"A": count_nt(bckgrnd, "G", "A"),
						"C": count_nt(bckgrnd, "G", "C"),
						"G": count_nt(bckgrnd, "G", "G"),
						"T": count_nt(bckgrnd, "G", "T"),	}

a_minus_anzahl["T"] = {	"A": count_nt(bckgrnd, "T", "A"),
						"C": count_nt(bckgrnd, "T", "C"),
						"G": count_nt(bckgrnd, "T", "G"),
						"T": count_nt(bckgrnd, "T", "T"),	}

# Summe aus Anzahl aller vorkommenden Nukleotide c-
sum_count_c_minus = 0

for x in a_minus_anzahl:
	for y in a_minus_anzahl[x]:
		sum_count_c_minus += a_minus_anzahl[x][y]

sum_count_c_minus = float(sum_count_c_minus)	#Typkonversion, um bei Berechungen genauere Zahlen zu erhalten


# Transitionsmatrix in Form von Dictionary a_minus berechnen
a_minus = {}

a_minus["A"] = {	"A": count_nt(bckgrnd, "A", "A") / sum_count_c_minus,
					"C": count_nt(bckgrnd, "A", "C") / sum_count_c_minus,
					"G": count_nt(bckgrnd, "A", "G") / sum_count_c_minus,
					"T": count_nt(bckgrnd, "A", "T") / sum_count_c_minus, }

a_minus["C"] = {	"A": count_nt(bckgrnd, "C", "A") / sum_count_c_minus,
					"C": count_nt(bckgrnd, "C", "C") / sum_count_c_minus,
					"G": count_nt(bckgrnd, "C", "G") / sum_count_c_minus,
					"T": count_nt(bckgrnd, "C", "T") / sum_count_c_minus, }

a_minus["G"] = {	"A": count_nt(bckgrnd, "G", "A") / sum_count_c_minus,
					"C": count_nt(bckgrnd, "G", "C") / sum_count_c_minus,
					"G": count_nt(bckgrnd, "G", "G") / sum_count_c_minus,
					"T": count_nt(bckgrnd, "G", "T") / sum_count_c_minus, }

a_minus["T"] = {	"A": count_nt(bckgrnd, "T", "A") / sum_count_c_minus,
					"C": count_nt(bckgrnd, "T", "C") / sum_count_c_minus,
					"G": count_nt(bckgrnd, "T", "G") / sum_count_c_minus,
					"T": count_nt(bckgrnd, "T", "T") / sum_count_c_minus, }

print("Transitionsmatrix a(-)-Modell")
print(a_minus)


#####################################
###Aufgabe 2: Mitochondriale Genom###
#####################################

sequence = list(fasta2array("test.fasta")[0])	#erste Array-element buchstabenweise gesplittet

	

#function zur Berechnung des likelihood-ratio s. fastafile: zu untersuchendes File; k: Startposition; w: Fensterbreite)  
def window_iteration(seq_list, k, w):
	res = []									#result list
	while (k < len(seq_list)-1-k+w-1):			#k muss kleiner sein als Laenge des Arrays minus Fensterbreite
		s = 0									# s: likelihood-ratio fuer entsprechendes Fenster
		i = k									# Variable i rueckt eine Position weiter, da k zuvor inkrementiert wurde
		while (i < k+w-1):
			s += log10(a_plus[seq_list[i]][seq_list[i+1]] / a_minus[seq_list[i]][seq_list[i+1]])
			i += 1
		res.append(s)
		k += 1
	return res

			
result_w51 = window_iteration(sequence,0,51)
result_w201 = window_iteration(sequence,0,201)
result_w401 = window_iteration(sequence,0,401)



######################################
###Aufgabe : Visualisierung im Plot###
######################################

list51= []
for x in range(len(result_w51)):
	list51.append(x)

plt.plot(list51, result_w51)
plt.show()


list201= []
for x in range(len(result_w201)):
	list201.append(x)

plt.plot(list201, result_w201)
plt.show()

list401= []
for x in range(len(result_w401)):
	list401.append(x)

plt.plot(list401, result_w401)
plt.show()





