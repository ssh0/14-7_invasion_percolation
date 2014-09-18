#! /usr/bin/env python
# -*- coding:utf-8 -*-
#
# written by Shotaro Fujimoto, August 2014.

from Tkinter import *
import sys
import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

class Percolation:

    def __init__(self, Lx=40, Ly=20):
        self.sub = None
        self.Lx = Lx
        self.Ly = Ly

    def perc_cluster(self):
        self.lattice = np.random.random([self.Lx+2, self.Ly])
        Lx = self.Lx
        Ly = self.Ly
        
        # 左端はすべて占有サイト
        self.lattice[:1,:] = 1.5
        
        self.lattice[Lx+1:,:] = 0
        if self.sub is None or not self.sub.winfo_exists():
            lattice = self.lattice
            ne = [(0, -1), (0, 1), (-1, 0), (1, 0)]
            
            # 周辺の点の座標を辞書形式で保持する
            nnsite = {(1, y): lattice[1, y] for y in range(Ly)}
            
            percolate = False
            while len(nnsite) != 0 and percolate == False:
                
                # 周辺の点で最も値の小さい格子の座標をmmに
                mm = min([(v, k) for k, v in nnsite.items()])[1]
                
                lattice[mm] = 1
                del nnsite[mm]
                
                # mmの周辺の点の座標をリストnnに(y方向に周期境界条件を適用)
                nn = [(mm[0] + nx, (mm[1] + ny)%Ly) for nx, ny in ne
                                if lattice[mm[0] + nx, (mm[1] + ny)%Ly] < 1]
                
                # nnの中で既にnnsiteに含まれているものを除く --> nnlist
                nnlist = list(set(nn) - set(nnsite.keys()))
                
                # nnsiteに新たな周辺の点を追加する
                for n in nnlist:
                    nnsite[n] = lattice[n]
                
                # 占有された格子点が右端である時，パーコレートしたと見なす
                if mm[0] == Lx:
                    percolate = True
            
            self.lattice = lattice[1:-1,:]
        return self.lattice
    
    def draw_canvas(self, rect, Lx, Ly):
        default_size = 640 # default size of canvas
        r = int(default_size/(2*Lx))
        fig_size_x = 2*r*Lx
        fig_size_y = 2*r*Ly
        margin = 10
        sub = Toplevel()
        
        sub.title('invasion percolation')
        self.canvas = Canvas(sub, width=fig_size_x+2*margin,
                    height=fig_size_y+2*margin)
        self.canvas.create_rectangle(margin, margin,
                    fig_size_x+margin, fig_size_y+margin,
                    outline='black', fill='white')
        self.canvas.pack()
        
        c = self.canvas.create_rectangle
        
        site = np.where(rect == 1)
        for m, n in zip(site[0], site[1]):
            c(2*m*r+margin, 2*n*r+margin,
                        2*(m+1)*r+margin, 2*(n+1)*r+margin,
                        outline='black', fill='black')
    
    def get_fractal_dim(self, trial=20, Lmin=20, Lmax=40, Lsample=10):
        
        # LminからLmaxの間の整数値で，できるだけlogにしたとき等間隔になるように
        L = np.array([int(i) for i 
                    in np.logspace(np.log10(Lmin), np.log10(Lmax), Lsample)])
        
        M_L = []
        for l in L:
            self.Lx = l*2
            self.Ly = l
            m_L = 0
            for i in range(trial):
                lattice = self.perc_cluster()
                
                # 中心のL×L格子中の占有サイト数を合計
                m_L += np.sum(lattice[int(l/2)+1:l+int(l/2),:] == 1)
            
            M_L.append(m_L/float(trial))
            print "L = %d, M_L = %f" % (l, M_L[-1])
        
        M_L = np.array(M_L)

        def fit_func(parameter0, L, M_L):
            log = np.log
            c1 = parameter0[0]
            c2 = parameter0[1]
            residual = log(M_L) - c1 - c2*log(L)
            return residual
        
        parameter0 = [0.1, 2.0]
        result = optimize.leastsq(fit_func, parameter0, args=(L, M_L))
        c1 = result[0][0]
        D = result[0][1]
        print "D = %f" % D
        
        def fitted(L, c1, D):
            return np.exp(c1)*(L**D)
        
        fig = plt.figure("Fractal Dimesion")
        ax = fig.add_subplot(111)
        ax.plot(L, M_L, '-o', label=r"$M(L)$")
        ax.plot(L, fitted(L, c1, D), label="fit func: D = %f" % D)
        ax.set_xlabel(r'$\ln L$', fontsize=16)
        ax.set_ylabel(r'$\ln M(L)$', fontsize=16)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ymargin(0.05)
        fig.tight_layout()
        plt.legend(loc='best')
        plt.show()

class TopWindow:
    
    def quit(self):
        self.root.destroy()
        sys.exit()
        
    def show_window(self, title="title", *args):
        self.root = Tk()
        self.root.title(title)
        frames = []
        for i, arg in enumerate(args):
            frames.append(Frame(self.root, padx=5, pady=5))
            for k, v in arg:
                Button(frames[i], text=k, command=v).pack(expand=YES, fill='x')
            frames[i].pack(fill='x')
        f = Frame(self.root, padx=5, pady=5)
        Button(f, text='quit', command=self.quit).pack(expand=YES, fill='x')
        f.pack(fill='x')
        self.root.mainloop()

if __name__ == '__main__':
    Lx = 40
    Ly = 20
    top = TopWindow()
    per = Percolation(Lx, Ly)
    count = 1

    def pr():
        global count
        d = per.canvas.postscript(file="figure_%d.eps" % count)
        print "saved the figure to a eps file"
        count += 1

    def pushed():
        per.perc_cluster()
        per.draw_canvas(per.lattice, Lx, Ly)

    def b4_pushed():
        trial = 100
        Lmin = 20
        Lmax = 100
        Lsample = 10
        per.get_fractal_dim(trial, Lmin, Lmax, Lsample)
    
    run = (('run', pushed), ('save canvas to sample.eps', pr))
    run2 = (('calculate D', b4_pushed),)
    top.show_window("Invasion Percolation", run, run2)

