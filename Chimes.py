__author__ = "Paolo Di Domenico"
__copyright__ = "Copyright 2017, Paolo Di Domenico"
__credits__ = ["Paolo Di Domenico"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Paolo Di Domenico"
__email__ = "yakkoroma@gmail.com"
__status__ = "Draft"

import math

class Atom:
    
    subshell_label = ['s','p','d','f','g','h','i','j','k']
    types = ['non metal','alkali metal','alkali earth metal','earth metal','halogen',
             'noble gas','transition metal','post transition metal','metalloids',
             'lanthanoids','actinoids','not classified']

    def __init__(self, proton_number):
        # teorical labels
        
        self.subshell_order = [
            {'1s':0},
            {'2s':0},
            {'2p':0},
            {'3s':0},
            {'3p':0},
            {'4s':0},
            {'3d':0},
            {'4p':0},
            {'5s':0},
            {'4d':0},
            {'5p':0},
            {'6s':0},
            {'4f':0},
            {'5d':0},
            {'6p':0},
            {'7s':0},
            {'5f':0},
            {'6d':0},
            {'7p':0}
        ]
       
        self.proton_number = proton_number
        """
        @TODO add isotopes
        """
        self.electron_number = proton_number
        self.max_quantic_number = int(math.ceil(math.sqrt(float(self.electron_number)/2.0)))
        aufbaum = self.aufbaum()
        self.electronic_configuration = aufbaum['configuration']
        self.block = aufbaum['block']
        self.period = self.get_period()
        self.valence_electron = self.get_valence_electron()
        self.group = self.valence_electron
        self.type = self.get_type()
        self.valence = self.get_valence()
    
    def get_type(self):
        """
        Cascade get type of atom
        @TODO non metal and metalloids in group p
              lack a group classification.
              check if is true
        """      
        if self.proton_number in [1,6,7,8,15,16,34]:
            return  0
        elif self.proton_number in [5,14,32,33,51,52,84]:
             return  8
        elif self.block == 's' and self.group == 1:
            return 1
        elif self.block == 's' and self.group == 2:
            return  2
        elif self.block == 'p' and self.group == 3:
            return  3
        elif self.block == 'p' and self.group == 7:
            return 4
        elif self.block == 'p' and self.group == 8:
            return 5
        elif self.block == 'p' and self.group != 8:
            return 7
        elif self.block == 'd' or self.block == 'f':
            return 6
        elif self.proton_numer <= 71 and self.proton_number >= 56:
            return 9
        elif self.proton_numer <= 103 and self.proton_number >= 89:
            return  10
        else:
            return 11
     
    def get_valence(self):
        """
        Get the valence as a list
        with the ecception for:
            B (5), Al(13);Si (14);O(8);F(9)
        """
        if self.group == 1:
            return [1]
        elif self.group == 2:
            return [2]
        elif self.group == 3:
            if self.proton_number in [5, 13]:
                return [3]
            else:
                return [1,3]
        elif self.group == 4:
            if self.proton_number == 14:
                return [4]
            else: 
                return [2,4]
        elif self.group == 5:
            return [3,5]
        elif self.group == 6:
            if self.proton_number == 8:
                return [2]
            else:
                return [2,4,6]
        elif self.group == 7:
            if self.proton_number == 9:
                return [1]
            else:
                return [1,3,5,7]
        else:
            """
            Noble gas
            """
            return 0

    def get_period(self):
        """
        The period is the maximum quantic number
        This function is added to give an idea
        that the period is the last occupied level
        """
        period = 0
        for i in self.electronic_configuration:
            # i.keys()[0] is a string, the first element is the main quantic number
            if int(i.keys()[0][0]) > period:
                period = int(i.keys()[0][0])
        return period

    def get_valence_electron(self):
        """
        The number of valence's electron is the sum
        of the electron in the last energetic level
        - which is the period -
        """
        valence = 0
        for i in self.electronic_configuration:
            if int(i.keys()[0][0]) == self.period:
                valence += i[i.keys()[0]]
        return valence

    def aufbaum(self):
        # computation coast is here is O(max_quantic_number ^ 2)
        n = self.max_quantic_number
        o, electron_sum = 0, 0
        slot_available = []
        """
        Build a list of dictionary of the avialable slot in orbital
        for the atom
        o is the main quantic number 
        l is the second quantic number (the type of orbital)
        """
        while o <= n:
            l = 0
            while l < o:
                """
                max_len is the number of electron the orbitals in 
                o energy level (2, one  for each spin) 
                Pauli exclusion rule
                Not contempling Hund rule (spin quantic number)
                an electron occupy before a free orbital and then the one 
                with the opposite spin
                """
                max_len = 2 * (2*l +1)          
                slot_available.append({str(o) + self.subshell_label[l]: max_len})
                l += 1
            o += 1
        """
        Order the slot for energetic level, fixed in a list
        """        
        for i in slot_available:
            for j, el in enumerate(self.subshell_order):
                if el.keys()[0] == i.keys()[0]:
                    self.subshell_order[j][i.keys()[0]] = i[i.keys()[0]]        
        slot_available = self.subshell_order       
        slot_taken = []
        sum_of_electron = 0
        """
        Calculate the configuration
        """
        for i, num in enumerate(slot_available):
            electron_number_in_slot = num[num.keys()[0]]
            if sum_of_electron <= self.electron_number:
                teorical_sum = sum_of_electron + electron_number_in_slot               
                if teorical_sum <= self.electron_number:
                    sum_of_electron += electron_number_in_slot
                    slot_taken.append(num)
                else:                                       
                    real_sum = self.electron_number - sum_of_electron
                    num[num.keys()[0]] = real_sum
                    if real_sum > 0:
                        slot_taken.append(num)
                        sum_of_electron += real_sum
        """
        Clean not taken slot
        """
        index = len(slot_taken)
        for i, l in enumerate(slot_taken):            
            if l[l.keys()[0]] == 0:
                index = i
                break        
        """
        Return the used orbital
        The last one name give the block
        """
        configuration = slot_taken[:index]        
        # i.keys()[0] is a string, the second element is the secondary quantic number
        return {'configuration':configuration, 'block': configuration[-1].keys()[0][1]}
        

    def show(self):
        print ("The atomic number is %s (no isotpes)" % (self.proton_number))
        print ("Its electronic configuration is %s" % self.electronic_configuration)
        print ("Its period is %s and the block is %s" % (self.period, self.block))
        print ("It has %s valence electron" % self.valence_electron)
        print ("It has valence: %s" % self.valence)
        print ("Its type is %s" % self.types[self.type])
