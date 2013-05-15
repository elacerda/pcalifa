'''
Created on 14/05/2013

@author: lacerda
'''

import pystarlight.io
import numpy as np
import atpy
from scipy import linalg
import PCAlifa

class ssp:
    ''' classe usada para ler e manipular os dados da base Starlight, 
        para ler os dados do CALIFA basta modificar esta classe que a
        classe GALAXY nao deve ser alterada.
        __init__() eh chamado toda a vez que vc instancia a classe:
        base = ssp(tabela)
    '''
    def __init__(self, table):
        self.nZ = len(np.unique(table.Z_base))
        self.nages = len(np.unique(table.age_base))
        self.shape = self.nZ, self.nages, -1
        self.l = table.l_ssp[0]
        self.f__Ztl = table.f_ssp.reshape([self.nZ, self.nages, -1])
        self.age_base = table.age_base.reshape([self.nZ, self.nages])[0]
        self.Z_base = table.Z_base.reshape([self.nZ, self.nages])[:, 0]

        self.Z_mask = (self.Z_base == 0.02)
        self.l_mask = (self.l > 3800) & (self.l < 6850)
        self.young_mask = (np.log10(self.age_base) < 7)
        self.old_mask = (np.log10(self.age_base) > 9.5)

        self.l_ssp = self.l[self.l_mask]

        self.pop_young()
        self.pop_old()
        self.set_lnorm(5635.)

    def pop_young_norm(self):
        self.f_young_norm = self.f_young / self.f_young[self.lnorm_i]

    def pop_old_norm(self):
        self.f_old_norm = self.f_old / self.f_old[self.lnorm_i]

    def pop_young(self):
        self.f_young = self.f__Ztl[self.Z_mask, self.young_mask][-1, self.l_mask]

    def pop_old(self):
        self.f_old = self.f__Ztl[self.Z_mask, self.old_mask][0, self.l_mask]

    def set_lnorm(self, l):
        self.lnorm_mask = (self.l_ssp >= l)
        self.lnorm = self.l_ssp[self.lnorm_mask][0]
        self.lnorm_i = np.where(self.l_ssp >= 1)[0][0]

        self.pop_young_norm()
        self.pop_old_norm()

class galaxy:
    ''' 
    coordenadas imagem: 
    y = linha
    x = coluna
    i = indice vetor imagem 
    imagem armazenada como um vetor com SIDE * SIDE posicoes 
    '''
    def __init__(self, side, base, rint, rpop):
        self.init = True
        self.side = side

        ''' criando a base com os espectros base '''
        self.base = ssp(base)
        self.fyoung__l = self.base.f_young
        self.fold__l = self.base.f_old
        self.fyoung_norm__l = self.base.f_young_norm
        self.fold_norm__l = self.base.f_old_norm
        self.nl = len(self.fyoung__l)


        ''' garantindo que o lado tenha tamanho impar '''
        if not (self.side % 2):
            self.side = side + 1
            print 'galaxy: new cube side: %d' % self.side

        ''' num de elementos '''
        self.n = np.int(self.side ** 2.)

        ''' acerta as coordenadas do nucleo para o centro do vetor'''
        self.set_nucleus_coords(self.n / 2)

        ''' calcula o raio de cada posicao xy'''
        self.set_radius()

        ''' calcula o vetor populacao(raio) usando como variavel controle RPOP '''
        self.set_rpop(rpop)

        ''' calcula o vetor de intensidades(raio) usando como variavel controle RINT '''
        self.set_rint(rint)

        ''' Cria o espectro de modo que:
            f__il[i, :] = int[i] * (x_y[i] * f_y + x_o[i] * f_o)
            f_norm__il[i, :] = int[i] * (x_y[i] * f_norm_y + x_o[i] * f_norm_o)
        '''
        self.set_f()

        self.init = False

    def set_nucleus_coords(self, i):
        self.i_nuc = i
        self.nuc_y, self.nuc_x = self.get_coords(self.i_nuc)

        if not self.init:
            self.set_radius()
            self.set_rpop(self.rpop)

    def set_radius(self, norm = False):
        self.r_norm = norm
        self.r = np.zeros(self.n, dtype = np.float64)

        for i in range(len(self.r)):
            self.r[i] = self.get_radius(i)

        if norm:
            if norm == 'max':
                self.r = self.r / self.r.max()
            else:
                self.r = self.r / self.r.sum()

    def set_rpop(self, rpop):
        self.rpop = rpop
        self.X()

        if not self.init:
            self.set_f()

    def set_rint(self, rint):
        self.rint = rint
        self.intens()

        if not self.init:
            self.set_f()

    def intens(self):
        self.intensity = (1. / np.exp(self.r / self.rint))

    def set_f(self):
        dim = np.ones((self.n, self.nl))
        self.f__il = ((self.Xyoung * dim.T).T * self.fyoung__l) + ((self.Xold * dim.T).T * self.fold__l)
        self.f__il *= (self.intensity * dim.T).T
        self.f_norm__il = ((self.Xyoung * dim.T).T * self.fyoung_norm__l) + ((self.Xold * dim.T).T * self.fold_norm__l)
        self.f_norm__il *= (self.intensity * dim.T).T

    def X(self):
        self.Xyoung = (1. / np.exp(self.r / self.rpop))
        self.Xold = 1. - self.Xyoung

    def get_coords(self, i):
        y = i / self.side
        x = i % self.side + 1

        return y, x

    def get_radius(self, i):
        y, x = self.get_coords(i)

        return ((1. * (x - self.nuc_x)) ** 2.0 + (1. * (y - self.nuc_y)) ** 2.0) ** 0.5


if __name__ == '__main__':
    basedir = '/home/lacerda/BC03models/'
    PadovaChab1994 = atpy.Table(basefile = basedir + 'Base.bc03.Padova1994.chab.All',
                                basedir = basedir,
                                type = 'starlightv4_base',
                                read_basedir = True)

    gal = galaxy(side = 51, base = PadovaChab1994, rint = 5., rpop = 5.)

    ''' Do what u want with GAL ;-) '''

    ''' PCA with GAL :P '''
    I__il, ms__l, covMat__ll, eigVal__k, eigVec__lk = PCAlifa.PCA(gal.f__il, gal.nl, 0)
    tomo__ik = np.dot(I__il, eigVec__lk)

    I_norm__il, ms_norm__l, covMat_norm__ll, eigVal_norm__k, eigVec_norm__lk = PCAlifa.PCA(gal.f_norm__il, gal.nl, 0)
    tomo__ik = np.dot(I_norm__il, eigVec_norm__lk)


