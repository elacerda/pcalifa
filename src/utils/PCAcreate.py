'''
Created on 14/05/2013

@author: lacerda
'''
import PCAlifa as pca
import pystarlight.io
import numpy as np
import atpy
from scipy import linalg
from matplotlib import pyplot as plt

def plotRebuildSpec(l, O, R, eVal, eVec, eVecUsed, npref):
    ''' criando uma string com os eigenvectors usados para reconstruir o cubo'''
    eVecUsedStr = ''
    ne = len(eVecUsed)

    for i in range(ne):
        eVecUsedStr = eVecUsedStr + '%d ' % eVecUsed[i]

    eVecUsedStr = eVecUsedStr.replace(eVecUsedStr[-2:], eVecUsedStr[-2])

    res = O - R
    adev = 100. * (1. / len(l)) * (np.abs(res) / O).sum()

    sigmaNReb = 0.
    sigmaReb = np.sqrt(eVal[:ne].sum())

    if ((ne + 1) != len(eVec)):
        sigmaNReb = np.sqrt(eVal[(ne + 1):].sum())

    sigmaRatio = sigmaNReb / sigmaReb

    f = plt.figure(figsize = (19.8, 10.8))
    ax = f.gca()
    ax.plot(l, O, label = 'O')
    ax.plot(l, R, label = 'R')
    ax.plot(l, res, label = 'res')
    textStr = 'adev =  %.4f %% - sigmaReb = %.2e - sigmaNReb = %.2e - ratio = %.4f' % (adev, sigmaReb, sigmaNReb, sigmaRatio)

    ax.text(0.01, 0.92, textStr, fontsize = 10, transform = ax.transAxes,
            horizontalalignment = 'left', verticalalignment = 'center', multialignment = 'left')

    plt.title('Rebuild using eigenvector(s): %s' % eVecUsedStr)
    ax.legend()
    ax.grid()
    f.savefig('%s_eVec-%s_rebSpec.png' % (npref, eVecUsedStr.replace(' ', '-')))
    plt.close()

class ssp:
    ''' classe usada para ler e manipular os dados da base Starlight, 
        para ler os dados do CALIFA basta modificar esta classe que a
        classe GALAXY nao deve ser alterada.
        __init__() eh chamado toda a vez que vc instancia a classe:
        base = ssp(tabela)
    '''
    def __init__(self, table):
        self.init = True
        self.nZ = len(np.unique(table.Z_base))
        self.nages = len(np.unique(table.age_base))
        self.l = table.l_ssp[0]
        self.f__Ztl = table.f_ssp.reshape([self.nZ, self.nages, -1])
        self.age_base = table.age_base.reshape([self.nZ, self.nages])[0]
        self.Z_base = table.Z_base.reshape([self.nZ, self.nages])[:, 0]

        ''' mascaras '''
        self.Z_mask = (self.Z_base >= 0.019) & (self.Z_base < 0.021)
        self.l_mask = (self.l > 3800) & (self.l < 6850)
        self.young_mask = (np.log10(self.age_base) < 7)
        self.old_mask = (np.log10(self.age_base) > 9.5)

        self.l_ssp = self.l[self.l_mask]
        self.set_lnorm(5635.)
        self.pop_young()
        self.pop_old()
        self.pop_young_norm()
        self.pop_old_norm()
        self.init = False

    def pop_young(self):
        self.f_young = self.f__Ztl[self.Z_mask, self.young_mask][-1, self.l_mask]

    def pop_old(self):
        self.f_old = self.f__Ztl[self.Z_mask, self.old_mask][0, self.l_mask]

    def pop_young_norm(self):
        self.f_young_norm = self.f_young / self.f_young[self.lnorm_mask][0]

    def pop_old_norm(self):
        self.f_old_norm = self.f_old / self.f_old[self.lnorm_mask][0]

    def set_lnorm(self, l):
        self.lnorm_mask = (self.l_ssp >= l)
        self.lnorm = self.l_ssp[self.lnorm_mask][0]

        if not self.init:
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
        self.fyoung_norm__l = self.base.f_young_norm
        self.fold_norm__l = self.base.f_old_norm
        self.l = self.base.l_ssp
        self.nl = len(self.l)

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
        '''
        XYOUNG tem dimensao de i
        XYOUNG * DIM.Transpose'''
        dim = np.ones((self.n, self.nl))

        self.f__il = ((self.Xyoung * dim.T).T * self.fyoung_norm__l) + ((self.Xold * dim.T).T * self.fold_norm__l)
        self.f__il *= (self.intensity * dim.T).T

    def X(self):
        self.Xold = (1. / np.exp(self.r / self.rpop))
        self.Xyoung = 1. - self.Xold

    def get_coords(self, i):
        y = i / self.side
        x = i % self.side + 1

        return y, x

    def get_radius(self, i):
        y, x = self.get_coords(i)

        return ((1. * (x - self.nuc_x)) ** 2.0 + (1. * (y - self.nuc_y)) ** 2.0) ** 0.5

if __name__ == '__main__':
    basedir = '/home/lacerda/BC03models/'
    PadovaSalp2000 = atpy.Table(basefile = basedir + 'Base.bc03.Padova2000.salp.All',
                                basedir = basedir,
                                type = 'starlightv4_base',
                                read_basedir = True)

    gal = galaxy(side = 51, base = PadovaSalp2000, rint = 5., rpop = 5.)

    ''' NOW Do what u want with GAL ;-) '''

##########################################################
    ''' PCA with GAL :P '''

    P = pca.PCAlifa()

    I__il, ms__l, covMat__ll, eigVal__k, eigVec__lk = P.PCA(gal.f__il, gal.nl, 0)
    ''' TOMO_KYX nao eh calculado pois nao existe P.K.ZoneToYX.'''
    tomo__ik, tomo__kyx = P.tomogram(I__il, eigVec__lk)
    ''' caso necessario podemos calcular TOMO__KYX usando GAL.SIDE '''
    tomo__kyx = tomo__ik.reshape(gal.side, gal.side, -1).T

    ''' plotando os 5 primeiros tomogramas e autoespectros '''
    P.tomoPlot(0, gal.l, tomo__kyx, eigVec__lk, eigVal__k, 'GalPadSalp2000_')
    P.tomoPlot(1, gal.l, tomo__kyx, eigVec__lk, eigVal__k, 'GalPadSalp2000_')
    P.tomoPlot(2, gal.l, tomo__kyx, eigVec__lk, eigVal__k, 'GalPadSalp2000_')
    P.tomoPlot(3, gal.l, tomo__kyx, eigVec__lk, eigVal__k, 'GalPadSalp2000_')
    P.tomoPlot(4, gal.l, tomo__kyx, eigVec__lk, eigVal__k, 'GalPadSalp2000_')

    ''' correlacoes '''
    by_radius = np.argsort(gal.r)

    colArr = [
            gal.Xyoung[by_radius],
            gal.Xold[by_radius],
            gal.intensity[by_radius],
    ]

    nRows = 5
    nCols = len(colArr) + 1
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.corrPlot(colArr[j], tomo__ik[:, i][by_radius], axArr[i, j])

        axArr[i, nCols - 1].plot(gal.l, eigVec__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    axArr[0, 0].set_title('Xyoung')
    axArr[0, 1].set_title('Xold')
    axArr[0, 2].set_title('Intensity')

    plt.suptitle('Correlations PC0 ... PC4 - OBS NORM')
    f.savefig('GalPadSalp2000_corre_obs_norm_0-4.png')
    plt.close()

    ''' reconstruindo cubos '''
    mask = np.zeros_like(eigVal__k, dtype = np.bool)
    mask[0] = True
    eigVecUsed = np.where([x for x in mask])[0]
    eigVecUsedStr = ''
    ne = len(eigVecUsed)

    for i in range(ne):
        eigVecUsedStr = eigVecUsedStr + '%d ' % eigVecUsed[i]

    eigVecUsedStr = eigVecUsedStr.replace(eigVecUsedStr[-2:], eigVecUsedStr[-2])

    ''' reconstrucao '''
    I_reconstr__il, f_reconstr__il = P.rebuildSpectra(tomo__ik[:, mask], eigVec__lk[:, mask], ms__l)

    ''' plot de espectro de algumas zonas '''
    for i in [0, gal.i_nuc, gal.i_nuc + (gal.side / 4)]:
        npref = 'GalPadSalp2000_zone-%04d' % i

        f = plt.figure(figsize = (19.8, 10.8))
        ax = f.gca()

        P.plotAxisZoneRebuildSpec(ax, gal.l, gal.f__il[i, :], f_reconstr__il[i, :],
                                  eigVal__k, eigVec__lk, mask, npref, resid = False)

        f.savefig('%s_eVec-%s_rebSpec.png' % (npref, eigVecUsedStr.replace(' ', '-')))
        plt.close()

    ''' reconstruindo usando outros eigenvectors '''
    mask = np.zeros_like(eigVal__k, dtype = np.bool)
    mask[1] = True
    eigVecUsed = np.where([x for x in mask])[0]
    eigVecUsedStr = ''
    ne = len(eigVecUsed)

    for i in range(ne):
        eigVecUsedStr = eigVecUsedStr + '%d ' % eigVecUsed[i]

    eigVecUsedStr = eigVecUsedStr.replace(eigVecUsedStr[-2:], eigVecUsedStr[-2])

    ''' reconstrucao '''
    I_reconstr__il, f_reconstr__il = P.rebuildSpectra(tomo__ik[:, mask], eigVec__lk[:, mask], ms__l)

    for i in [0, gal.i_nuc, gal.i_nuc + (gal.side / 4)]:
        npref = 'GalPadSalp2000_zone-%04d' % i

        f = plt.figure(figsize = (19.8, 10.8))
        ax = f.gca()

        P.plotAxisZoneRebuildSpec(ax, gal.l, gal.f__il[i, :], f_reconstr__il[i, :],
                                  eigVal__k, eigVec__lk, mask, npref, resid = False)

        ax.set_ylabel('PC %s' % eigVecUsedStr)
        f.savefig('%s_eVec-%s_rebSpec.png' % (npref, eigVecUsedStr.replace(' ', '-')))
        plt.close()

    ''' reconstruindo usando eVec 0, 1 '''
    mask = np.zeros_like(eigVal__k, dtype = np.bool)
    mask[[0, 1]] = True
    eigVecUsed = np.where([x for x in mask])[0]
    eigVecUsedStr = ''
    ne = len(eigVecUsed)

    for i in range(ne):
        eigVecUsedStr = eigVecUsedStr + '%d ' % eigVecUsed[i]

    eigVecUsedStr = eigVecUsedStr.replace(eigVecUsedStr[-2:], eigVecUsedStr[-2])

    ''' reconstrucao '''
    I_reconstr__il, f_reconstr__il = P.rebuildSpectra(tomo__ik[:, mask], eigVec__lk[:, mask], ms__l)

    for i in [0, gal.i_nuc, gal.i_nuc + (gal.side / 4)]:
        npref = 'GalPadSalp2000_zone-%04d' % i

        f = plt.figure(figsize = (19.8, 10.8))
        ax = f.gca()

        P.plotAxisZoneRebuildSpec(ax, gal.l, gal.f__il[i, :], f_reconstr__il[i, :],
                                  eigVal__k, eigVec__lk, mask, npref, resid = False)

        f.savefig('%s_eVec-%s_rebSpec.png' % (npref, eigVecUsedStr.replace(' ', '-')))
        plt.close()
