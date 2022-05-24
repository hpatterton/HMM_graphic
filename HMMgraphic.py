import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad
from collections import Counter

class HMM:

    def __init__(self, number_nodes, molecule_type):
        self.dictionary_of_transitions = {'M0_M1':0,'M0_I0':1,'M0_D1':2,'I0_M1':3,'I0_I0':4,'I0_D1':5,'MX_MY':6,'MX_IX':7,'MX_DY':8,'IX_MY':9,'IX_IX':10,'IX_DY':11,'DX_MY':12,'DX_IX':13,'DX_DY':14,'ML_Mend':15,'ML_IL':16,'IL_Mend':17,'IL_IL':18,'DL_Mend':19,'DL_IL':20}
        self.dictionary_of_emmission_nodes = {'M':0,'I':1}
        self.dictionary_of_amino_acids = {'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
        self.dictionary_of_nucleotides = {'G':0,'A':1,'T':2,'C':3}
        # 0 is nucleotide/aa, 1 is node number (0 is start, N is end), 2 is M or I
        self.molecule_type = molecule_type
        if self.molecule_type == 'P': # make sure it is not a 'shallow' list
            self.array_of_emmisions = [[[0 for i in range(2)] for j in range(number_nodes+2)] for k in range(len(self.dictionary_of_amino_acids))]
        else:
            self.array_of_emmisions = [[[0 for i in range(2)] for j in range(number_nodes+2)] for k in range(len(self.dictionary_of_nucleotides))]
        print(self.array_of_emmisions)
        self.number_of_nodes = number_nodes
        self.array_of_transitions = [[0]*21]*(number_nodes+2) # 21 is the sum of all possible transitions listed in self.dictionary_of_transitions
        self.array_of_nodes = [[0],[0],[0]]
        self.node_type_indexes = {'M':0,'I':1,'D':2}
        self.make_end_node()
        for node_number in range(number_nodes):
            for node_type in range(len(self.node_type_indexes)):
                self.make_node([], node_type)
        self.make_start_node()
        return

    def set_transition_value(self,value, node_index, sub_node_index):
        self.array_of_transitions[node_index][sub_node_index] = value
        return

    def get_transition_value(self, node_index, sub_node_index):
        return(self.array_of_transitions[node_index][sub_node_index])

    def make_end_node(self):
        node = HMM_node([],'M',[-1,-1,-1])
        self.array_of_nodes[self.node_type_indexes['M']][0] = node
        #self.array_of_nodes[0][self.node_type_indexes['I']] = 0
        return(node)

    def make_node(self, emmission_array, node_type_index):
        # number of M nodes currently
        index_of_this_node = len(self.array_of_nodes[0])
        # # point to earlier nodes
        if node_type_index == self.node_type_indexes['I']:
            child_pointers = [index_of_this_node-1,-1,index_of_this_node-1]
        elif node_type_index == self.node_type_indexes['M'] or node_type_index == self.node_type_indexes['D']:
            child_pointers = [index_of_this_node-1,index_of_this_node-1,index_of_this_node-1]
        node = HMM_node([], node_type_index, child_pointers)
        self.array_of_nodes[node_type_index].append(node)
        return

    def make_start_node(self):
        index_of_this_node = len(self.array_of_nodes[0])
        child_pointers = [index_of_this_node-1, index_of_this_node-1, index_of_this_node-1]
        node = HMM_node([], self.node_type_indexes['M'], child_pointers)
        self.array_of_nodes[self.node_type_indexes['M']].append(node)
        return

    def list_nodes(self):
        #print(self.array_of_nodes)
        number_of_nodes = self.array_of_nodes[0]
        print(number_of_nodes)
        # print(self.array_of_nodes[0][0].get_node_info())
        # print(self.array_of_nodes[1][0].get_node_info())
        # print(self.array_of_nodes[1][1].get_node_info())
        # print(self.array_of_nodes[1][2].get_node_info())
        # for node in range(0,number_of_nodes):
        #     for i in range(0,3):
        #         if self.array_of_nodes[node][i] > 0:
        #             print(node, self.array_of_nodes[node][i].get_node_info())
        return

    def Read_FastA(path):
        f = open(path, 'r')
        name = []
        sequence = []
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                name.append(line[1:].rstrip('\n'))
                if len(current_sequence) > 0:
                    sequence.append(current_sequence)
                current_sequence = ''
            else:
                current_sequence += line.rstrip('\n')
        sequence.append(current_sequence)
        return (name, sequence)







class HMM_node:

    def __init__(self, list_of_symbols, type_of_node, indexes_of_child_nodes):
        #self.emmission_symbols = list_of_symbols
        self.node_type = type_of_node
        self.child_indexes = indexes_of_child_nodes
        return

    def get_node_type(self):
        print(self.node_type)
        return

    def get_node_info(self):
        print('type of mode = ',self.node_type, "child pointers=",self.child_indexes)
        return

class HMM_sequence:

    def __init__(self):
        self.number_of_sequences = 0
        self.start_block_is_conserved = True
        self.name = []
        self.sequence = []
        self.sequence_array = np.array([['A','-','-','-','-'],['A','-','A','-','A'],['G','A','A','A','C'],['G','A','T','A','C'],['G','T','G','A','C']])
        self.sequence_blocks = [[],[]]
        self.order_of_blocks = []
        self.conservation_cutoff = 0.75
        self.dictionary_of_amino_acids = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
                                          'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16,
                                          'V': 17, 'W': 18, 'Y': 19}
        self.dictionary_of_nucleotides = {'G': 0, 'A': 1, 'T': 2, 'C': 3}
        return

    def Read_FastA(self, path):
        try:
            with open(path, 'r') as f:
                pass
        except FileNotFoundError:
            print('File not found')
            return False
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                self.name.append(line[1:].rstrip('\n'))
                if len(current_sequence) > 0:
                    self.sequence.append(list(current_sequence))
                current_sequence = ''
            else:
                current_sequence += line.rstrip('\n')
        self.sequence.append(list(current_sequence))
        self.sequence_array = np.array(self.sequence)
        if len(self.sequence) > 0:
            number_of_sequences = len(self.sequence)
        else:
            return False
        sequence_length = self.sequence[0]
        result = True
        for i in range(1,number_of_sequences):
            if sequence_length == len(self.sequence[i]):
                result = False
        if result == False:
            return False

        return

    def get_pseudocounts(self, numerator, denominator, molecule_type):
        if molecule_type == 'P':
            return((numerator+1)/(denominator+20))
        else:
            return((numerator+1)/(denominator+4))

    def get_emmission_list(self,list_of_symbols, molecule_type):
        if molecule_type == 'P':
            emmision_array = [0 for i in range(len(self.dictionary_of_amino_acids))]
            symbol_count = Counter(list_of_symbols)
            total_symbols = sum(symbol_count.values())
            for character in self.dictionary_of_amino_acids:
                numerator = symbol_count[character]
                denominator = total_symbols
                # You will need to call pseudocounts here if you want to use it
                emmision_array[self.dictionary_of_amino_acids[character]] = numerator/denominator
        elif molecule_type == 'D':
            emmision_array = emmision_array = [0 for i in range(len(self.dictionary_of_nucleotides))]
            symbol_count = Counter(list_of_symbols)
            total_symbols = sum(symbol_count.values())
            for character in self.dictionary_of_nucleotides:
                numerator = symbol_count[character]
                denominator = total_symbols
                # You will need to call pseudocounts here if you want to use it
                entry = self.dictionary_of_nucleotides[character]
                emmision_array[entry] = numerator#/denominator
        else:
            emmision_array = []

        return emmision_array

    



    def is_conserved(self, row):
        row_as_list = row.tolist()
        length_of_block = len(row_as_list)
        number_of_hyphens = row_as_list.count('-')
        if (length_of_block-number_of_hyphens)/length_of_block >= self.conservation_cutoff:
            return True
        else:
            return False


    def identify_blocks(self):
        start = 0
        current_column = start
        length = len(self.sequence_array[0,:].tolist())
        while current_column < length:
            while (current_column < length) and (self.is_conserved(self.sequence_array[:,current_column])):
                current_column += 1
            if (start < length) and (current_column <= length) and start != current_column:
                self.sequence_blocks[0].append((start,current_column))
                #print(start,current_column)
                self.order_of_blocks.append(0)
                start = current_column
            while (current_column < length) and not(self.is_conserved(self.sequence_array[:,current_column])):
                current_column += 1
            if (start < length) and (current_column <= length):
                self.sequence_blocks[1].append((start,current_column))
                #print(start, current_column)
                self.order_of_blocks.append(1)
                start = current_column
        #print(self.order_of_blocks)

        return

    def number_of_conserved_nodes(self):
        sum=0
        length = len(self.sequence_blocks[0])
        for i in range(length):
            sum += self.sequence_blocks[0][i][1]-self.sequence_blocks[0][i][0]
        return(sum)

class HMM_draw:

    def __init__(self):
        self.radius=1
        self.scale=0.25
        self.text_size = 20*2*self.radius*self.scale
        self.x_origin = 2
        self.y_origin = 2
        self.vertical_spacing = 4
        self.horizontal_spacing = 4
        self.number_of_nodes = 10
        self.figure_size_x = 16
        self.figure_size_y = 9
        self.fig = plt.figure(figsize=(16, 9))
        self.ax = self.fig.add_axes((0, 0, 1, 1))
        self.ax.set_xlim(0, 16)
        self.ax.set_ylim(0, 9)
        #ax.set_facecolor(bg_color)

        self.ax.tick_params(bottom=False, top=False,
                       left=False, right=False)
        self.ax.tick_params(labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
        return

    def draw_circle(self,x,y):
        theta = np.linspace(0, 2 * np.pi, 100)
        print(np.sin(theta))
        self.ax.plot(self.scale*(x*self.horizontal_spacing + self.radius * np.cos(theta))+self.x_origin, self.scale*(y*self.vertical_spacing + self.radius * np.sin(theta))+self.y_origin,color="midnightblue",)
        return

    def draw_square(self,x,y):
        self.ax.plot((self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin,self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin),(self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin,self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin),color="midnightblue",linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin), color="midnightblue" ,lw=1,linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin), color="midnightblue" ,linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin), color="midnightblue" ,linestyle='-')
        return

    def draw_diamond(self,x,y):
        self.ax.plot((self.scale*(x*self.horizontal_spacing-self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing)+self.x_origin), (self.scale*(y*self.vertical_spacing)+self.y_origin, self.scale*(y*self.vertical_spacing+self.radius)+self.y_origin), color="midnightblue",
                     linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing)+self.x_origin, self.scale*(x*self.horizontal_spacing+self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing+self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing)+self.y_origin), color="midnightblue", lw=1,
                     linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing+self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing)+self.x_origin), (self.scale*(y*self.vertical_spacing)+self.y_origin, self.scale*(y*self.vertical_spacing-self.radius)+self.y_origin), color="midnightblue",
                     linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing)+self.x_origin, self.scale*(x*self.horizontal_spacing-self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing-self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing)+self.y_origin), color="midnightblue",
                     linestyle='-')
        self.draw_text(self.x_origin + self.scale * (x * self.horizontal_spacing),
                       self.y_origin + self.scale * (y * self.vertical_spacing), 'I' + str(x), 'center', 'center', self.text_size,
                       'black')

        return

    def show_drawing(self):
        plt.show()
        return

    def draw_begin(self):
        # draw 'Begin'
        x=0
        y=0
        self.draw_square(x, y)
        self.draw_main_match_block_arrows(x,y)
        self.draw_text(self.x_origin, self.y_origin, 'Begin', 'center', 'center', self.text_size, 'black')

        # draw insert diamond
        self.draw_diamond(x, y+1)
        self.draw_main_block_diamond_arrows(x, y+1)
        return

    def draw_end(self):
        self.draw_square((self.number_of_nodes+1), 0)
        self.draw_text(self.x_origin + self.scale * ((self.number_of_nodes+1) * self.horizontal_spacing),
                       self.y_origin + self.scale * (0 * self.vertical_spacing), 'End', 'center',
                       'center', self.text_size,
                       'black')
        return

    def draw_text(self,x,y,text, h_align, v_align, font_size, col):
        text_kwargs = dict(ha=h_align, va=v_align, fontsize=font_size, fontname='Arial', color=col)
        plt.text(x, y, text, **text_kwargs)
        return

    def draw_arrow(self,x1,y1,x2,y2):
        self.ax.annotate("",(x1,y1),(x2,y2),arrowprops=dict(arrowstyle="<|-"))
        return

    def draw_main_block_diamond_arrows(self,x,y):
        self.draw_insert_block_arrows(self.scale*(x*self.horizontal_spacing-self.radius)+self.x_origin, self.scale*(y*self.horizontal_spacing)+self.y_origin)
        #draw to next square
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing - 0.5 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius), self.y_origin + self.scale * (0 * self.vertical_spacing+self.radius))
        # draw to next delete circle
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing + 0.5 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius), self.scale*(2*self.vertical_spacing-0.5*self.radius)+self.y_origin)
        return

    def draw_last_block_diamond_arrows(self, x, y):
        self.draw_insert_block_arrows(self.scale * (x * self.horizontal_spacing - self.radius) + self.x_origin,
                                      self.scale * (y * self.horizontal_spacing) + self.y_origin)
        # draw to next square
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing - 0.5 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius), self.y_origin + self.scale * (0 * self.vertical_spacing+self.radius))
        return

    def draw_main_block_delete_arrows(self,x,y):
        self.draw_arrow(self.x_origin + self.scale*(x*self.horizontal_spacing+ 0.707 * self.radius), self.y_origin+self.scale*(y*self.vertical_spacing-0.707*self.radius),
                        self.x_origin + self.scale * ((x+1) * self.horizontal_spacing-0.5*self.radius), self.scale * (0 * self.vertical_spacing+self.radius)+self.y_origin)
        self.draw_arrow(self.x_origin + self.scale*(x*self.horizontal_spacing), self.y_origin+self.scale*(y*self.vertical_spacing-self.radius),
                        self.x_origin + self.scale*x*self.horizontal_spacing, self.y_origin+self.scale*((y-1)*self.vertical_spacing+self.radius))
        self.draw_arrow(self.x_origin + self.scale*(x*self.horizontal_spacing+self.radius), self.y_origin+self.scale*(y*self.vertical_spacing),
                        self.x_origin + self.scale*((x+1)*self.horizontal_spacing-self.radius), self.y_origin+self.scale*(y*self.vertical_spacing))
        return

    def draw_end_block_delete_arrows(self, x, y):
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.707 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing - 0.707 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing-0.5*self.radius), self.scale * (0 * self.vertical_spacing+self.radius)+self.y_origin)
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing),
                        self.y_origin + self.scale * (y * self.vertical_spacing - self.radius),
                        self.x_origin + self.scale * x * self.horizontal_spacing,
                        self.y_origin + self.scale * ((y - 1) * self.vertical_spacing + self.radius))
        return

    def draw_main_match_block_arrows(self,x,y):
        self.draw_arrow(self.x_origin + self.scale * (x+self.radius), self.y_origin,
                        self.x_origin + self.scale * (x+self.horizontal_spacing - self.radius), self.y_origin)
        self.draw_arrow(self.x_origin + self.scale * x, self.y_origin + self.scale * self.radius,
                        self.x_origin + self.scale * x,
                        self.y_origin + self.scale * (self.vertical_spacing - self.radius))
        self.draw_arrow(self.x_origin+self.scale*(x+self.radius), self.y_origin + self.scale * self.radius,
                        self.x_origin + self.scale * (x+self.horizontal_spacing - 0.707 * self.radius),
                        self.y_origin + self.scale * (2 * self.vertical_spacing - 0.707 * self.radius))
        return

    def draw_last_main_match_block_arrows(self, x, y):
        self.draw_arrow(self.x_origin + self.scale * (x+self.radius), self.y_origin,
                        self.x_origin + self.scale * (x+self.horizontal_spacing - self.radius), self.y_origin)
        self.draw_arrow(self.x_origin + self.scale * x, self.y_origin + self.scale * self.radius,
                        self.x_origin + self.scale * x,
                        self.y_origin + self.scale * (self.vertical_spacing - self.radius))
        return

    def draw_insert_block_arrows(self,x,y):
        self.draw_circular_arrow(x, y)
        return

    def draw_circular_arrow(self,x,y):
        #fig = plt.figure(figsize=(9, 9))
        #ax = plt.gca()
        x_radius = self.scale
        y_radius=self.scale
        base_angle = 0 # the 'reference' placement from where angles are measured is 0
        start_angle = 55 # the angle will start at 225 degrees, from a 45 degree angle
        end_angle =247 # the length of the arrow is 270 degrees from the start_angle
        arrow_capstyle = 'round'
        arrow_linestyle = '-'
        arrow_linewidth = 1
        arrow_colour = 'black'
        # ========Line
        arc = Arc([x, y], x_radius, y_radius, angle=start_angle, theta1=base_angle, theta2=end_angle, capstyle=arrow_capstyle, linestyle=arrow_linestyle, lw=arrow_linewidth, color=arrow_colour)
        self.ax.add_patch(arc)
        # ========Create the arrow head
        arrow_head_x = x + (x_radius / 2) * np.cos(rad(end_angle + start_angle))  # Do trig to determine end position
        arrow_head_y = y + (y_radius / 2) * np.sin(rad(end_angle + start_angle))
        # Create triangle as arrow head # (x,y) # number of vertices # radius # orientation
        self.ax.add_patch(RegularPolygon((arrow_head_x, arrow_head_y), 3, x_radius/15, rad(start_angle+end_angle), color=arrow_colour))
        return

    def draw(self):
        self.draw_begin()
        for i in range(self.number_of_nodes):
            self.draw_square((i+1), 0)
            self.draw_text(self.x_origin + self.scale * ((i+1) * self.horizontal_spacing),
                               self.y_origin + self.scale * (0 * self.vertical_spacing), 'M' + str(i+1), 'center',
                               'center', self.text_size,
                               'black')
            if i < self.number_of_nodes-1:
                self.draw_main_match_block_arrows((i+1)*self.horizontal_spacing, 0)
            else:
                self.draw_last_main_match_block_arrows((i+1)*self.horizontal_spacing, 0)
            self.draw_diamond((i + 1), 1)
            if i < self.number_of_nodes-1:
                self.draw_main_block_diamond_arrows((i+1), 1)
            else:
                self.draw_last_block_diamond_arrows(i+1, 1)
            self.draw_circle((i + 1), 2)
            self.draw_text(self.x_origin + self.scale * ((i + 1) * self.horizontal_spacing),
                           self.y_origin + self.scale * (2 * self.vertical_spacing), 'D' + str(i + 1), 'center',
                           'center', self.text_size,
                           'black')
            if i < self.number_of_nodes - 1:
                self.draw_main_block_delete_arrows(i+1, 2)
            else:
                self.draw_end_block_delete_arrows(i+1, 2)
            # if i < self.number_of_nodes-1:
            #      self.draw_main_block_delete_arrows(i, 2)
            # else:
            #     self.draw_end_block_delete_arrows(i, 2)
        self.draw_end()
        return

    def draw_figure_title(self, text):
        self.draw_text((self.scale*(self.number_of_nodes+2)*self.horizontal_spacing+self.x_origin)*0.5, self.scale*3 * self.vertical_spacing+self.y_origin, 'My HMM 1', 'center', 'center', self.text_size, 'black')











#HMM = HMM(3,'N')
#HMM.list_nodes()
# hmm_sequences = HMM_sequences()
# print(hmm_sequences.sequence_array)
# hmm_sequences.identify_blocks()
# print(hmm_sequences.sequence_blocks)
# print(hmm_sequences.number_of_conserved_nodes())
# print(HMM.number_of_nodes)
# print(HMM.array_of_nodes)
#HMM.make_end_node()
# hmm_draw = HMM_draw()
# hmm_draw.draw()
# hmm_draw.draw_figure_title('My HMM')
# #hmm_draw.draw_arrow(1,1,4,4)
# hmm_draw.show_drawing()
hmm_sequence = HMM_sequence()
result = hmm_sequence.get_emmission_list(['G','G','A','T','C','C','C',], 'D')
print(result)