# -*- coding: utf-8 -*-
# file: phase_diagram.py

# This code is part of blendpy.
# MIT License
#
# Copyright (c) 2025 Leandro Seixas Rocha <leandro.fisica@gmail.com> 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
import pandas as pd
from tqdm import tqdm
from ase.parallel import parprint as print


class PhaseDiagram():
    def __init__(self, enthalpy: np.ndarray, entropy: np.ndarray, temperatures: np.ndarray = np.arange(300, 3001, 50)):
        """
        Initialize the phase diagram with enthalpy, entropy, and optional temperature values.

        Parameters:
        enthalpy (np.ndarray): Array of enthalpy values.
        entropy (np.ndarray): Array of entropy values.
        temperatures (np.ndarray, optional): Array of temperature values. Defaults to np.arange(300, 3001, 50).

        Raises:
        NotImplementedError: If the lengths of enthalpy and entropy arrays do not match.
        """
        if np.issubdtype(enthalpy.dtype, np.number) == False:
            raise TypeError("Enthalpy should be a numeric array.")
        if np.issubdtype(entropy.dtype, np.number) == False:
            raise TypeError("Entropy should be a numeric array.")
        if np.issubdtype(temperatures.dtype, np.number) == False:
            raise TypeError("Temperatures should be a numeric array.")
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.temperatures = temperatures
        self.gibbs_free_energy = None
        self.npoints = len(enthalpy)
        self.x = np.linspace(0,1,self.npoints)


    def get_gibbs_free_energy(self) -> np.ndarray:
        """
        Calculate and return the Gibbs free energy.

        The Gibbs free energy is computed using:
        G = H - T * S
        where:
        - G is the Gibbs free energy
        - H is the enthalpy
        - T is the temperature
        - S is the entropy

        Returns:
            np.ndarray: The calculated Gibbs free energy.
        """
        self.gibbs_free_energy = self.enthalpy - self.temperatures[:, None] * self.entropy
        return self.gibbs_free_energy


    def get_spinodal_decomposition(self, verbose: bool = True) -> pd.DataFrame:
        if verbose:
            print("-----------------------------------------------")
            print("  Spinodal curve")
            print("-----------------------------------------------")

        spinodal = []
        for t in self.temperatures:
            dx = 1/(self.npoints-1)
            diff2_gibbs = np.gradient(np.gradient(self.gibbs_free_energy, dx), dx)
            idx = np.argwhere(np.diff(np.sign(diff2_gibbs))).flatten()
            data = [t, self.x[idx]]
            flattened_array = np.concatenate([np.atleast_1d(item) for item in data])
            spinodal.append(flattened_array)
        
        # df_spinodal = (
        #     pd.concat([
        #         pd.DataFrame(spinodal)[[1, 0]].rename(columns={1: "x", 0: "t"}),
        #         pd.DataFrame(spinodal)[[2, 0]].rename(columns={2: "x", 0: "t"}).iloc[::-1].reset_index(drop=True)
        #     ], ignore_index=True)
        #     .dropna()
        # )
        # df0 = pd.DataFrame(spinodal, columns=["t", "a", "b"])
        # df1 = df0[["a","t"]]
        # df1.columns = ["x","t"]
        # df2 = df0[["b","t"]]
        # df2.columns = ["x","t"]
        # reversed_df2 = df2.iloc[::-1].reset_index(drop=True)
        # df_result = pd.concat([df1, reversed_df2], axis=0, ignore_index=True)
        # df_spinodal = df_result.dropna()
        
        # # spinodal critical point
        # (x_c, T_c) = df_spinodal.iloc[df_spinodal['t'].argmax()]
        # print("    Spinodal critical point (x_c, T_c):", (x_c, T_c))
        
        return spinodal
        # return df_spinodal


    def get_new_spinodal_decomposition(self, verbose: bool = True) -> pd.DataFrame:

        # Calculate the first derivative with respect to composition
        dG_dcomp = np.gradient(self.gibbs_free_energy, self.x, axis=1)

        # Calculate the second derivative with respect to composition
        d2G_dcomp2 = np.gradient(dG_dcomp, self.x, axis=1)

        # Find the compositions where the second derivative is zero or changes sign.
        spinodal_points = []

        for i, T in enumerate(self.temperatures):
            row = d2G_dcomp2[i, :]
            # Loop over adjacent composition points
            for j in range(len(self.x) - 1):
                # Check if exactly zero (or nearly zero)
                if np.isclose(row[j], 0, atol=1e-6):
                    spinodal_points.append((T, self.x[j]))
                # Check for a sign change between row[j] and row[j+1]
                elif row[j] * row[j+1] < 0:
                    # Linear interpolation to estimate the zero crossing composition:
                    c1, c2 = self.x[j], self.x[j+1]
                    d1, d2 = row[j], row[j+1]
                    # Linear interpolation formula: c_zero = c1 - d1*(c2-c1)/(d2-d1)
                    c_zero = c1 - d1 * (c2 - c1) / (d2 - d1)
                    spinodal_points.append((c_zero, T))

        # Convert the list of tuples into a pandas DataFrame
        spinodal_df = pd.DataFrame(spinodal_points, columns=['x', 't'])
        return spinodal_df

    
    def _dif_gibbs(self):
        """
        Calculate the gradient of the Gibbs free energy with respect to composition.

        This method computes the numerical gradient of the Gibbs free energy
        with respect to composition using a central difference method.

        Returns:
            np.ndarray: A 2D array containing the gradient of the Gibbs free energy
            with respect to composition. The shape of the array is (n_temperatures, npoints).
        """
        dx = 1 / (self.npoints - 1)
        dif_gibbs = np.gradient(self.gibbs_free_energy, dx, axis=1)
        return dif_gibbs


    def _gibbs_taylor(self, x0: int):
        """
        Compute the Taylor expansion of the Gibbs free energy around a given point x0.
        Parameters:
        -----------
        x0 : int
            The index of the point around which to perform the Taylor expansion.
        Returns:
        --------
        taylor : numpy.ndarray
            A 2D array of shape (n_temperatures, npoints) containing the Taylor expansion
            of the Gibbs free energy around x0 for each temperature.
        """
        x = self.x  # 1D array with shape (npoints,)
        g = self.gibbs_free_energy  # 2D array with shape (n_temperatures, npoints)
        dg = self._dif_gibbs()  # 2D array with shape (n_temperatures, npoints)
        
        # For each temperature (each row), compute the Taylor expansion around x0:
        # g[:, x0] is an array of shape (n_temperatures,)
        # We add a new axis to make it (n_temperatures, 1) so that it broadcasts with x - x[x0] (which is (npoints,))
        taylor = g[:, x0][:, None] + dg[:, x0][:, None] * (x - x[x0])
        return taylor  # Shape: (n_temperatures, npoints)


    def _is_convex(self, x0: int):
        g = self.gibbs_free_energy()    # shape: (n_temperatures, npoints)
        taylor = self._gibbs_taylor(x0)   # shape: (n_temperatures, npoints)

        # Check if g >= taylor for every x in each temperature (row)
        convexity_at_t = np.all(g >= taylor, axis=1)
        return convexity_at_t  # returns a boolean array of shape (n_temperatures,)
    

    def _miscibility_gap(self):
        convexity_flags = []
        for x0 in self.x:
            for ti in self.temperatures:
                for xi in self.x:
                    convexity_flags.append(self._is_convex(x0)[ti,xi])
        return convexity_flags
    
        # # Ensure endpoints are considered convex.
        # convexity_flags[0] = convexity_flags[-1] = True

        # if not all(convexity_flags):
        #     first_idx = convexity_flags.index(False)
        #     last_idx = len(convexity_flags) - 1 - convexity_flags[::-1].index(False)
        #     return (x[first_idx], x[last_idx])
        # return (None, None)


    # def get_binodal_curve(self) -> pd.DataFrame:
    #     print("-----------------------------------------------")
    #     print("\033[36mBinodal curve\033[0m")
    #     print("-----------------------------------------------")

        
    #     # Compute the miscibility gap for each temperature.
    #     binodal_data = [self._miscibility_gap(T, A, B, slope, eps, npoints) for T in temperatures]
    #     df_gap = pd.DataFrame(binodal_data, columns=['xi', 'xf'])
    #     df = pd.DataFrame({'t': temperatures})
    #     df = pd.concat([df, df_gap], axis=1).dropna().reset_index(drop=True)
        
    #     # Prepare lower and upper halves of the solvus curve.
    #     df_lower = df[['xi', 't']].copy()
    #     df_lower.columns = ["x", "t"]
    #     df_upper = df[['xf', 't']].copy()
    #     df_upper.columns = ["x", "t"]
    #     df_upper = df_upper.iloc[::-1].reset_index(drop=True)
        
    #     # Concatenate to form the complete solvus curve.
    #     df_binodal = pd.concat([df_lower, df_upper], ignore_index=True)

    #     (x_c, T_c) = df_binodal.iloc[df_binodal['t'].argmax()]
    #     print("    Binodal critical point (x_c, T_c):", (x_c, T_c))
        
    #     return df_binodal
