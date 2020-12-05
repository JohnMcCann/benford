#!/usr/bin/env python3

import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
cc = plt.rcParams['axes.prop_cycle'].by_key()['color']

from utils.radix import *


class benford:
    _ordinal_indicator = {1:'st', 2:'nd', 3:'rd'}
    def __init__(self, data=None, frequency=None, n_digits=1, fd_loc=1, base=10):
        """
        Description:
            Object for analysising data with respects to Benford's law of
            anomalous numbers.
            
        Keyword arguments:
            data: data to be analysized, raw numbers
            frequency: frequency of digits of the raw numbers (counts)
            n_digits: number of digits
            df_loc: location of first digit of n_digit chunk (can be < 0)
            base: base/radix being used
        """
        # Error check arguments
        if base < 2 or type(base) != int:
            print('Error: Base must be integer larger than 1.')
            return
        self.n_digits = n_digits
        self.fd_loc = fd_loc
        self.base = base
        # Determine numbers of interest (leading digits in radix 10)
        self.numbers_10 = np.asarray(
            range(self.base**(self.fd_loc+self.n_digits-2),
                  self.base**(self.fd_loc+self.n_digits-1), 1), dtype=int
        )
        # Determine numbers of interest in radix-b
        upb = int(self.base**self.n_digits)
        if self.fd_loc == 1:
            upb -= 1
        self.numbers_b = [baseB_to_str(base10_to_baseB(n, self.base)[0]
                                       [-self.n_digits:], self.base)
                          for n in self.numbers_10[:upb]]
        # Calculate probability distribution
        self.calc_pdf()
        if data is not None:
            self.data = np.asarray(data)
            # Drop data that doesn't conform to request
            self.data = self.data[self.data >= self.numbers_10[0]]
            self.build_freq()
            self.ndata = self.frequency.sum()
        if frequency is not None:
            if ((len(frequency) != self.base**self.n_digits and self.fd_loc != 1)
                or (len(frequency) != self.base**self.n_digits-1 and self.fd_loc == 1)):
                print('ERROR: Frequency length inconsistent with base '
                      'and first digit location.')
                return
            self.frequency = frequency
            self.ndata = self.frequency.sum()
        return


    def calc_pdf(self):
        """
        Description:
            Calculations the pdf of Benford's Law based on the radix, number
            of digits considered, and location of first digit.
        """
        if self.fd_loc == 1:
            self.pdf = np.log(1.+1./self.numbers_10)/np.log(self.base)
        else:
            self.pdf = np.zeros(self.base**self.n_digits)
            for n in self.numbers_10:
                i = base10_to_baseB(n, self.base)[0]
                i = i[-self.n_digits:]
                i = baseB_to_base10(i, self.base)
                self.pdf[i] += np.log(1.+1./n)/np.log(self.base)
        return


    def build_freq(self):
        """
        Description:
            Takes raw data and calculates the frequency of the string of
            digits under consideration for Benford's Law.
        """
        if self.fd_loc == 1:
            self.frequency = np.zeros(self.base**(self.n_digits-1)
                                      *(self.base-1),
                                      dtype=int)
        else:
            self.frequency = np.zeros(self.base**self.n_digits, dtype=int)
        lob = self.fd_loc-1
        upb = lob+self.n_digits
        for n in self.data:
            i = baseB_to_base10(
                base10_to_baseB(n, self.base)[0][lob:upb], self.base)
            if self.fd_loc == 1:
                i -= self.base**(self.n_digits)
            self.frequency[i] += 1
        # Also set proportion from frequency
        self.proportion = self.frequency/len(self.data)
        return


    def goodness(self):
        N = len(self.data)
        # Person's Chi squared test
        self.chi2 = N*((self.proportion-self.pdf)**2/self.pdf).sum()
        self.chi2_alpha = 1.-chi2.cdf(self.chi2, len(self.frequency)-1)
        # Cramer-von Mises statistics
        self.S_CvM = self.proportion.cumsum()
        self.T_CvM = self.pdf.cumsum()
        self.Z_CvM = self.S_CvM - self.T_CvM
        self.t_CvM = np.asarray(list((self.pdf[:-1]+self.pdf[1:])/2)
                               +[(self.pdf[-1]+self.pdf[0])/2])
        self.Z_CvM_mean = (self.t_CvM*self.Z_CvM).sum()
        # Cramer-von Mises test
        self.W_CvM = N*(self.t_CvM*self.Z_CvM**2).sum()
        self.U_CvM = N*(self.t_CvM*(self.Z_CvM-self.Z_CvM_mean)**2).sum()
        self.A_CvM = N*(self.t_CvM*self.Z_CvM**2
                        /(self.T_CvM*(1-self.T_CvM)))[:-1].sum()
        return


    def sim_conf(self, alpha=0.05, method='Goodman'):
        """
        Description:
            Calculates the simultaneous confidence intervals (SCI) or
            confidence band (CB)
            
        Keyword arguments:
            alpha: significance level
            method: See reference for methods
            
        Reference:
            https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0151235 
            
        Returns:
            tuple of the (lower, upper) simultaneous confidence intervals
        """
        if self.frequency is None:
            print('ERROR: Must provide data or frequency before calculating'
                  ' error bars.')
            return
        N_data = self.ndata
        dof = len(self.frequency)
        if method == 'Goodman' or method == 'Quesenberry':
            if method == 'Quesenberry':
                c = chi2.ppf(1-alpha, dof)
            elif method == 'Goodman':
                c = chi2.ppf(1-alpha/dof, 1)
            knot, shift = (np.zeros(dof) for i in range(2))
            for i in range(dof):
                knot[i] = (c+2*self.frequency[i])/(2*(N_data+c))
                shift[i] = np.sqrt(
                    c*(c+4*self.frequency[i]*(N_data-self.frequency[i])/N_data)
                )/(2*(N_data+c))
        return (knot, shift)


    def null_test(self, alpha=0.05, method='Goodman', verbose=True):
        """
        Description:
            Given a confidence band determined by the chosen method, checks
            if the null hypothesis is rejected.
        
        Keyword arguments:
            alpha: significance level
            method: method used for sim_conf() call
            
        Returns:
            Array of which digits rejected null hypothesis.
        """
        knot, shift = self.sim_conf(alpha=alpha, method=method)
        reject = np.asarray([False,]*len(self.pdf))
        rejected = None
        for i, p in enumerate(self.pdf):
            if not(p <= knot[i]+shift[i] and p >= knot[i]-shift[i]):
                reject[i] = True
        if verbose:
            if reject.any():
                print(f'Null hypothesis rejected for alpha={alpha:g} using '
                      f"the {method:s} method.")
                if self.base == 10:
                    rejected = [n for n, c in zip(self.numbers_10, reject) if c]
                    print(f"Digits rejected: {rejected}")
                else:
                    rejected = [n for n, c in zip(self.numbers_b, reject) if c]
                    print(f"Digits rejected: {', '.join(rejected):s}")
            else:
                print(f'Null hypothesis cannot be rejected for alpha={alpha:g} '
                      f"using the {method:s} method.")
        return rejected


    def plot_benford_data(self, ax, normalize=True, alpha=0.05,
                          method='Goodman', print_numbers=False,
                          sparse_frame=False):
        """
        Description:
            Make boiler template plot for visualizing data in light of
            Benford's Law of anomalous numbers.
            
        Arguments:
            ax: matplotlib axis object to plot

        Keyword arguments:
            normalize: Plot data in terms of proportions (boolean)
            method: method used for calculating confidence band

        Returns:
            pdf line object, and bar and errorbar container object
        """
        # Plot distribution, Benford's law, and error bars
        xticks = range(1, len(self.frequency)+1, 1)
        knot, shift = self.sim_conf(alpha=alpha, method=method)
        if normalize:
            pdf = ax.plot(xticks, self.pdf, 'rx', ms=15, mew=2,
                          label="Benford's Law")
            bars = ax.bar(xticks, self.frequency/self.ndata, label='Sample')
            err = ax.errorbar(xticks, knot, yerr=shift, ecolor='m', ls='',
                              capsize=10, capthick=3, lw=3,
                              label=method+' 95% CB')
            ax.set_ylabel('Proportion')
        else:
            pdf = ax.plot(xticks, self.ndata*self.pdf, 'rx', ms=15, mew=2,
                          label="Benford's Law")
            bars = ax.bar(xticks, self.frequency, label='Sample')
            err = ax.errorbar(xticks, self.ndata*knot, yerr=self.ndata*shift,
                              ecolor='m', ls='', capsize=10, capthick=3, lw=3,
                              label=method+' 95% CB')
            ax.set_ylabel('Frequency')
        # Display numbers at bottom of bars
        if print_numbers:
            for bar in bars:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    0.01,
                    round(bar.get_height(), 2),
                    fontsize=15,
                    horizontalalignment='center',
                    color='white',
                    weight='bold'
                )
        if sparse_frame:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            # ax.spines['bottom'].set_color('#DDDDDD')
            ax.tick_params(bottom=False, left=False)
            ax.set_axisbelow(True)
            ax.yaxis.grid(True)#, color='#EEEEEE')
            ax.xaxis.grid(False)
        # Set xtick labels
        if self.base == 10:
            ax.set_xticklabels(self.numbers_b)
        else:
            ax.set_xticklabels(self.numbers_b, rotation=90)
        ax.legend()
        ax.set_xticks(xticks)
        ax.set_xlabel(f'Digit')
        title = "Benford's Law:\n"
        n = self.fd_loc
        title += rf"${n:d}^{{{self._ordinal_indicator.get(n,'th')}}}$ "
        if self.n_digits != 1:
            n += self.n_digits-1
            title += (
                rf"thru ${n:d}^{{{self._ordinal_indicator.get(n,'th')}}}$ "
            )
        title += rf"digit in radix-{self.base:d} (N={self.ndata})"
        ax.set_title(title)
        return (pdf, bars, err)


    def confidence_plot(self, ax, alphas=[0.0001, 0.001, 0.01, 0.05]):
        alphas.sort()
        xticks = range(1, len(self.frequency)+1, 1)
        l1, = ax.plot(xticks, [0,]*len(xticks), 'o', color='w', ms=5,
                      markerfacecolor='none')
        l2, = ax.plot(xticks, self.proportion-self.pdf, '+', color='k', ms=15,
                      mew=3)
        # Legend for zeroes and sample data points
        legend = ax.legend([l1, l2], ['Zeroes', 'Sample'])
        ax.add_artist(legend)
        lines, labels = ([] for i in range(2))
        for i, alpha in enumerate(alphas):
            knot, shift = self.sim_conf(alpha=alpha)
            label = f"{100*(1-alpha):g}% CB"
            l = ax.fill_between(xticks, knot-shift-self.pdf, knot+shift-self.pdf,
                                 alpha=1, color=cc[i])
            lines.append(l)
            labels.append(label)
        ax.set_xlabel('Digit')
        ax.set_ylabel('Difference')
        ax.set_title("Proportion minus Benford's Law")
        ax.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -.225),
                  fancybox=False, shadow=False, ncol=min(4,len(alphas)))
        # Set xtick labels
        if self.base == 10:
            ax.set_xticklabels(self.numbers_b)
        else:
            ax.set_xticklabels(self.numbers_b, rotation=90)
        ax.set_xticks(xticks)
        return