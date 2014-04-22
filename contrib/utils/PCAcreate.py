'''
Created on 14/05/2013

@author: lacerda

TODO: 
    2013-05-30 - Doc das funcs.

'''
import PCAlifa as PCA
import pystarlight.io
import numpy as np
import atpy
from scipy import linalg
from matplotlib import pyplot as plt
from astropy.io import ascii

def plotRebuildSpec(l, O, R, eVal, eVec, eVecUsed, npref):
    ''' criando uma string com os eigenvectors usados para reconstruir o cubo '''
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
    ''' 
        Classe usada para ler e manipular os dados da base Starlight, 
        para ler os dados do CALIFA basta modificar esta classe que a
        classe GALAXY nao deve ser alterada.
    '''
    def __init__(self, table, lc = [3800, 6850]):
        self.init = True

        self.nZ = len(np.unique(table.Z_base))
        self.nages = len(np.unique(table.age_base))
        self.l = table.l_ssp[0]
        self.f__Ztl = table.f_ssp.reshape([self.nZ, self.nages, -1])
        self.age_base = table.age_base.reshape([self.nZ, self.nages])[0]
        self.Z_base = table.Z_base.reshape([self.nZ, self.nages])[:, 0]
        self.Z_mask = (self.Z_base >= 0.019) & (self.Z_base < 0.021)
        self.setLambdaConstrains(lc)
        self.setAgeMask(6.5, 9.5)
        self._setVars()

        self.init = False

    def _setVars(self):
        self.l_ssp = self.l[self.l_mask]
        self.setLambdaNorm(5635.)
        self.pop_young()
        self.pop_old()
        self.pop_young_norm()
        self.pop_old_norm()

    def setLambdaConstrains(self, lc = []):
        lc = np.array(lc)
        s = np.argsort(lc)
        ldown = lc[s][0]
        lup = lc[s][1]

        self.lc = lc[s]
        self.l_mask = (self.l >= ldown) & (self.l <= lup)

        if not self.init:
            self._setVars()

    def setAgeMask(self, ageYoung, ageOld):
        self.young_mask = (np.log10(self.age_base) < ageYoung)
        self.old_mask = (np.log10(self.age_base) > ageOld)

        if not self.init:
            self._setVars()

    def pop_young(self):
        self.f_young = self.f__Ztl[self.Z_mask, self.young_mask][-1, self.l_mask]

    def pop_old(self):
        self.f_old = self.f__Ztl[self.Z_mask, self.old_mask][0, self.l_mask]

    def pop_young_norm(self):
        self.f_young_norm = self.f_young / self.f_young[self.lnorm_mask][0]

    def pop_old_norm(self):
        self.f_old_norm = self.f_old / self.f_old[self.lnorm_mask][0]

    def setLambdaNorm(self, l):
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
            i = indice pixel 
        imagem armazenada como um vetor com SIDE * SIDE posicoes 
    '''
    def __init__(self, base, side = 1, rpop = 0., rint = 0., rdust = 0., Xo0 = 0., tauV0 = 0., extLawFile = False):
        self.init = True

#        self.rpop = rpop
#        self.rint = rint
#        self.rdus = rdust
#        self.Xo0 = Xo0
#        self.tauV0 = tauV0

        self.init_vars()

        self.set_side(side)

        if base:
            ''' criando a base com os espectros base '''
            lc = [3800, 6850]
            self.set_base(base, lc)

            ''' calcula o vetor populacao(raio) usando como variavel controle RPOP '''
            self.set_rpop(rpop, Xo0)

            ''' calcula o vetor de intensidades(raio) usando como variavel controle RINT '''
            self.set_rint(rint)

            ''' calcula o vetor de poeira(raio) usando como varival rdust e tauV0) '''
            self.set_rdust(rdust, tauV0, extLawFile)

            ''' calcula todos os espectros da galaxia  '''
            self.set_f()

        self.init = False

    def init_vars(self):
        self.extLaw = False
        self.dim__li = False
        self.XpopOld__i = False
        self.XpopOld__il = False
        self.XpopYoung__i = False
        self.XpopYoung__il = False
        self.fold_norm__l = False
        self.fyoung_norm__l = False
        self.intensity__i = False
        self.intensity__il = False
        self.tauV__i = False
        self.tauV__il = False
        self.q_CCM__l = False
        self.r__i = False

    def set_side(self, side):
        if side == 0:
            side = 1

        self.side = side

        ''' num de elementos '''
        self.n = np.int(self.side ** 2.)

        ''' acerta as coordenadas do nucleo para o centro do vetor'''
        self.set_nucleus_coords(self.n / 2)

        ''' garantindo que o lado tenha tamanho impar '''
        if not (self.side % 2):
            self.side = side + 1
            print 'galaxy: new cube side: %d' % self.side

        ''' calcula o raio de cada posicao xy'''
        self.set_radius()

    def set_nucleus_coords(self, i):
        self.i_nuc = i
        self.nuc_y, self.nuc_x = self.get_coords(self.i_nuc)

        if not self.init:
            self.set_radius()

            if (self.extLaw):
                self.set_rpop(self.rpop, self.Xo0)

    def set_radius(self, norm = False):
        self.r_norm = norm
        self.r__i = np.zeros(self.n, dtype = np.float64)

        for i in range(len(self.r__i)):
            self.r__i[i] = self.get_radius(i)

        if norm:
            if norm == 'max':
                self.r__i = self.r__i / self.r__i.max()
            else:
                self.r__i = self.r__i / self.r__i.sum()

    def set_base(self, base, lc = [3800, 6850]):
        self.base = ssp(base, lc)
        self.fyoung_norm__l = self.base.f_young_norm
        self.fold_norm__l = self.base.f_old_norm
        self.l = self.base.l_ssp
        self.nl = len(self.l)

        ''' matriz com 1 em todas as posicoes com dimensao nl x n '''
        self.dim__li = np.ones((self.nl, self.n))

        if not self.init:
            self.set_rpop(self.rpop, self.Xo0)
            self.set_rint(self.rint)
            self.set_rdust(self.rdust, self.tauV0)
            self.set_f()

    def set_rpop(self, rpop, Xo0):
        self.rpop = rpop
        self.Xo0 = Xo0
        self._Xpop()

        if not self.init:
            self.set_f()

    def _Xpop(self):
        self.XpopOld__i = self.Xo0 * np.exp(-1. * self.r__i / self.rpop)
        self.XpopOld__il = (self.XpopOld__i * self.dim__li).T
        self.XpopYoung__i = (1. - self.XpopOld__i)
        self.XpopYoung__il = (self.XpopYoung__i * self.dim__li).T

    def set_rint(self, rint):
        self.rint = rint
        self._intens()

        if not self.init:
            self.set_f()

    def _intens(self):
        self.intensity__i = np.exp(-1. * self.r__i / self.rint)
        self.intensity__il = (self.intensity__i * self.dim__li).T

    def set_rdust(self, rdust, tauV0, extLawFile = False):
        self.rdust = rdust
        self.tauV0 = tauV0

        if extLawFile:
            self.extLaw = self.readExtLaw(extLawFile)

        self._q()
        self._tau()

        if not self.init:
            self.set_f()

    def readExtLaw(self, file):
        self.extLawFile = file
        columns = ('lambda', 'q_CCM', 'q_CAL', 'q_HZ1', 'q_HZ2', 'q_HZ3', 'q_HZ4', 'q_HZ5', 'q_GD1', 'q_GD2', 'q_GD3')
        t = ascii.read(self.extLawFile,
                       Reader = ascii.FixedWidthNoHeader,
                       delimiter = ' ',
                       data_start = 0,
                       data_end = 19001,
                       names = columns)
        return t

    def _tau(self):
        self.tauV__i = self.tauV0 * np.exp(-1. * self.r__i / self.rdust)
        self.tauV__il = (self.tauV__i * self.dim__li).T

    def _q(self):
        if self.extLaw:
            l = np.array(self.extLaw['lambda'], dtype = np.float64)
            lmask = (l >= self.l[0]) & (l <= self.l[-1])
            self.q_CCM__l = np.array(self.extLaw['q_CCM'][lmask], dtype = np.float64)
        else:
            self.q_CCM__l = np.ones(self.nl, dtype = np.float64)

    def set_f(self):
        self.fint__il, self.fobs__il = self._f()

    def get_coords(self, i):
        y = i / self.side
        x = i % self.side + 1

        return y, x

    def get_radius(self, i):
        y, x = self.get_coords(i)

        return ((1. * (x - self.nuc_x)) ** 2.0 + (1. * (y - self.nuc_y)) ** 2.0) ** 0.5

    def set_no_dust(self):
        self.set_rdust(0., 0.)

    def _f(self):
        Xo = self.XpopOld__il
        fo = self.fold_norm__l
        Xy = self.XpopYoung__il
        fy = self.fyoung_norm__l
        Iint = self.intensity__il
        tauV = self.tauV__il
        q = self.q_CCM__l

        fint = Iint * ((Xo * fo) + (Xy * fy))
        fobs = fint * np.exp(-tauV * q)

        return fint, fobs

if __name__ == '__main__':
    basedir = '/home/lacerda/BC03models/'
    PadovaSalp2000 = atpy.Table(basefile = basedir + 'Base.bc03.Padova2000.salp.All',
                                basedir = basedir,
                                type = 'starlightv4_base',
                                read_basedir = True)

    gal = galaxy(base = PadovaSalp2000,
                 side = 21,
                 rpop = 10000.,
                 rint = 5.,
                 rdust = 1.,
                 Xo0 = 0.5,
                 tauV0 = 0.,
                 extLawFile = '/home/lacerda/CALIFA/ExtLaws.out',
                 )

    ''' NOW Do what u want with GAL ;-) '''

##########################################################
    ''' PCA with GAL :P '''
    P = PCA.PCAlifa()

    I__il, ms__l, covMat__ll, eigVal__k, eigVec__lk = P.PCA(gal.fobs__il, gal.nl, 0)
    ''' TOMO_KYX nao eh calculado pois nao existe P.K.ZoneToYX.'''
    tomo__ik, tomo__kyx = P.tomogram(I__il, eigVec__lk)
    ''' caso necessario podemos calcular TOMO__KYX usando GAL.SIDE '''
    tomo__kyx = tomo__ik.reshape(gal.side, gal.side, -1).T

    ''' plotando os 5 primeiros tomogramas e autoespectros '''
    for ti in range(5):
        P.tomoPlot(tomo__kyx, gal.l, eigVec__lk, eigVal__k, ti, 'GalPadSalp2000_')

    ''' correlacoes '''
    by_radius = np.argsort(gal.r__i)

    colArr = [
            gal.XpopYoung__i[by_radius],
            gal.XpopOld__i[by_radius],
            gal.intensity__i[by_radius],
            gal.tauV__i[by_radius]
    ]

    colNames = [
            r'$10^6\ to\ 10^{7.5}yr$',
            r'$10^{7.5}\ to\ 10^{8.5}yr$',
            r'$10^{8.5}\ to\ 10^{9.5}yr$',
            r'$10^{9.5}\ to\ 10^{10.2}yr$',
            r'eigenvector',
    ]

    nRows = 5
    nCols = len(colArr) + 1
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.correlationAxisPlot(colArr[j], tomo__ik[:, i][by_radius], axArr[i, j])

        axArr[i, nCols - 1].plot(gal.l, eigVec__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

    for i in range(nCols):
        axArr[0, i].set_title(colNames[i])

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

        P.zoneRebuildSpecAxisPlot(ax, gal.l, gal.fobs__il[i, :], f_reconstr__il[i, :], eigVal__k, eigVec__lk, mask, npref, 10, False)

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

        P.zoneRebuildSpecAxisPlot(ax, gal.l, gal.fobs__il[i, :], f_reconstr__il[i, :], eigVal__k, eigVec__lk, mask, npref, 10, False)

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

        P.zoneRebuildSpecAxisPlot(ax, gal.l, gal.fobs__il[i, :], f_reconstr__il[i, :], eigVal__k, eigVec__lk, mask, npref, 10, False)

        f.savefig('%s_eVec-%s_rebSpec.png' % (npref, eigVecUsedStr.replace(' ', '-')))
        plt.close()
