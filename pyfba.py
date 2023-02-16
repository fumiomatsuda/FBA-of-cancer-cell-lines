#!/usr/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        pyfba.py
# Purpose:     A simple module for flux balance analysis
#
# Author:      Fumio_Matsuda
#
# Created:     20/09/2022
# Copyright:   (c) Fumio_Matsuda 2022
# Licence:     MIT license
#-------------------------------------------------------------------------------

"""pyfba.py:MetabolicModel class

This module includes all class and its functions of pyfba.

The module includes::

    MetabolicModel class

Todo:
    * Write documents of all functions

"""
import re, csv
import pulp #conda install -c conda-forge pulp

class MetabolicModel:
    """Class of MetabolicModel

    An instances of this class has fuctions to execute FBA.

    Returns:
        instance: instance of metabolic model


    Examples:
        >>> model = pyfba.MetabolicModel(modelfilename, reversible = "<=>" , onedirection = "-->", plus = "+", output = "normal")

    Attributes (incomplete)::

            self.separator = {}
            self.constraint = {}
            self.optstep = {}
            self.inconstant = {}
            self.inconstant_array = []
            self.inconstant_onlyincrease = {}
            self.metabolites = {}
            self.reactions = {}
            self.reactionid_array = []
            self.maboliteid_array = []


    """

    def __init__(self, modelfilename, reversible = "<=>" , onedirection = "-->", plus = "+", output = "normal"):

        '''Generator of new instance
        Args:
            modelfilename (str): Path and name of metabolic model file

            reversible (str): Separator of reversivle reactions

            onedirection (str): Separator of irreversivle reactions

            plus (str): Separator of each metabolite in reactions

            mode (str): "normal" or "debug"

        '''


        format = "text"


        self.separator = {}
        self.constraint = {}
        self.optstep = {}
        self.inconstant = {}
        self.inconstant_array = []
        self.inconstant_onlyincrease = {}
        self.metabolites = {}
        self.reactions = {}
        self.reactionid_array = []
        self.maboliteid_array = []

        self.separator["reversible"] = reversible
        self.separator["onedirection"] = onedirection
        self.separator["plus"] = plus
        separator_onedirection = self.separator["reversible"]
        separator_reversible = self.separator["onedirection"]
        separator_plus = self.separator["plus"]


        #
        # Collection of metabolite names
        #
        reactionid = 0

        with open(modelfilename, 'r') as f:
            if format == "text":
                reader = csv.reader(f, delimiter='\t')
            elif format == "csv":
                reader = csv.reader(f, dialect='excel')
            else:
                print("Unknown format!")
                f.close()
                return False

            for i, row in enumerate(reader):
                if len(row) < 6:
                    continue

                rid = row[0].replace('"', "")		#ID
                if output == "debug":
                    print(rid)
                rid = self.check_reactionname(rid)


                reaction = row[5].replace('"', "") #	反応式
                flag = row[1].replace('"', "")
                group = row[2].replace('"', "") #反応の説明
                lb = row[3].replace('"', "") #下限
                ub = row[4].replace('"', "") #上限
                reaction = re.sub(r'^\s+',"",reaction)
                reaction = re.sub(r'\s+$',"",reaction)
                if flag == "0":
                    continue
                pat = re.compile(r'\s+{0}\s+|\s+{1}\s+'.format(separator_onedirection, separator_reversible))
                reactions = re.split(pat, reaction)
                if output == "debug":
                    print(rid, reaction, group, separator_onedirection, separator_reversible)
                if not len(reactions) == 2:
                    continue
                substrate = reactions[0]
                product = reactions[1]
                substrates = re.split(r'\s+\+\s+', substrate)
                products = re.split(r'\s+\+\s+', product)
                self.reactions[rid] = {}
                self.reactions[rid]["reaction"] = reaction
                self.reactions[rid]["group"] = group
                self.reactions[rid]["lb"] = lb
                self.reactions[rid]["ub"] = ub
                self.reactions[rid]["metabolites"] = {}
                self.reactions[rid]["reactionid"] = reactionid

                reactionid = reactionid + 1
                self.reactionid_array.append(rid)

                for compound in substrates:
                    cps =  re.split(r' +', compound)
                    #
                    # 1.1348 13BDglcn[c]
                    #
                    if re.match(r'^[0-9]+\.*[0-9]*$', cps[0]):
                        compound = cps[1]
                        compound = self.check_metabolitename(compound)
                        self.reactions[rid]["metabolites"][compound] = float(cps[0]) * -1.0
                    #
                    # (1.1348) 13BDglcn[c]
                    #
                    elif re.match(r'^\([0-9]+\.*[0-9]*\)$', cps[0]):
                        compound = cps[1]
                        compound = self.check_metabolitename(compound)
                        temp = str(cps[0])
                        temp = temp.replace('(','')
                        temp = temp.replace(')','')
                        self.reactions[rid]["metabolites"][compound] = float(temp) * -1.0
                    else:
                        compound = self.check_metabolitename(compound)
                        self.reactions[rid]["metabolites"][compound] = -1.0
                    self.metabolites[compound] = {"metaboliteid": 1, "reactions":{}}

                for compound in products:
                    cps =  re.split(r' +', compound)
                    if re.match(r'^[0-9]+\.*[0-9]*$', cps[0]):
                        compound = cps[1]
                        compound = self.check_metabolitename(compound)
                        self.reactions[rid]["metabolites"][compound] = float(cps[0]) * 1.0
                    #
                    # (1.1348) 13BDglcn[c]
                    #
                    elif re.match(r'^\([0-9]+\.*[0-9]*\)$', cps[0]):
                        compound = cps[1]
                        compound = self.check_metabolitename(compound)
                        temp = str(cps[0])
                        temp = temp.replace('(','')
                        temp = temp.replace(')','')
                        self.reactions[rid]["metabolites"][compound] = float(temp) * 1.0
                    else:
                        compound = self.check_metabolitename(compound)
                        self.reactions[rid]["metabolites"][compound] = 1.0
                    self.metabolites[compound] = {"metaboliteid": 1, "reactions":{}}


    def check_metabolitename(self, metid):

        """Method to check metabolie ID. Current version does nothing.
        13BDglcn[c] => 13BDglcn[c]

        Args:
            metid (str): Metabolite ID

        Returns:
            metid (str): Modified ID of metabolite

        Examples:
            >>> compound = self.check_metabolitename(compound)


        History:

        """

        if re.match(r'^[0-9]+', metid):
            metid = ""+metid
        return metid

    def check_reactionname(self, rid):
        """Method to check reaction ID.
        Letter "R" is added to any reaction IDs beginning with numbers.
            1RXM => R1RSN
        "-", "[", and "] are replaced with "_"
            1-glc[c] => R1_glc_c_

        Args:
            rid (str): Reaction ID

        Returns:
            rid (str): Modified Reaction ID

        Examples:
            >>> rid = self.check_reactionname(rid)


        History:

        """

        if re.match(r'^[0-9]+', rid):
            rid = "R"+rid
        rid = rid.replace('-', "_")		#ID
        rid = rid.replace('[', "_")		#ID
        rid = rid.replace(']', "_")		#ID
        return rid

    def solve(self, buffer = 0.9999, output = "normal"):
        """Method to perform liner optimization.
        Some optimization function is at least required to conduct.
        This method performs two step optimizations
        1st: Optimization of the objective funtion to find an optival value Oopt.
        2nd: Minimization of a net absolute metabolite flux level (Sigma(i)(abs(vi)) when the level of objective funtion is constrained at Oopt * buffer.
        "buffer" is a buffer to enable the 2nd minimization.


        See:
            add_optstep()

        Args:
            buffer (float): Reduce to 0.9999, 0.999, 0.99 when 2nd optimization is infeasible.
            output (str): "normal" or "debug". Detailed process is shown in the "debug" mode.

        Returns:
            status: "optimal" indicates that 2nd optimization was successfully finished
            objective (float): Optimization value of the optimization function

        Examples:
            >>> status,objective = model.solve()
        See


        History:

        """
        # Metabolite names
        self.maboliteid_array = list(sorted(self.metabolites.keys()))
        for i, met in enumerate(self.maboliteid_array):
            self.metabolites[met]["metaboliteid"] = i
        for i, reaction in enumerate(self.reactionid_array):
            for metabolite in self.reactions[reaction]["metabolites"]:
                self.metabolites[metabolite]["reactions"][reaction] =  self.reactions[reaction]["metabolites"][metabolite]
        #
        # PuLP LpProblem
        #
        self.problem = pulp.LpProblem('sample', pulp.LpMaximize)
        self.x = []
        x = []
        #
        # Reactions
        #
        for rid in self.reactionid_array:
            lb= self.reactions[rid]["lb"]
            ub= self.reactions[rid]["ub"]
            if output == "debug":
                print(rid, lb, ub)
            self.x.append(pulp.LpVariable(rid, float(lb), float(ub), pulp.LpContinuous))
        #
        # Objective function
        #
        coeffs = []
        rids = []
        for optstep in self.optstep:
            reactionid = self.reactions[optstep]["reactionid"]
            coef =  self.optstep[optstep]
            coeffs.append(float(coef))
            rids.append(self.x[reactionid])
            if output == "debug":
                print(coeffs, rids, reactionid)

        self.problem += pulp.lpDot(coeffs, rids)
        #self.problem += pulp.lpSum(coeffs * rids)
        #
        # Constrains, fixed value
        #
        for rid in self.constraint:
            reactionid = self.reactions[rid]["reactionid"]
            value = self.constraint[rid]
            if output == "debug":
                print(rid, reactionid, value)
            self.problem +=  self.x[reactionid] == float(value)

        #
        # Constrains, stoichiometry
        #
        for metabolite in self.maboliteid_array:
            coeffs = []
            rids = []
            for rid in self.metabolites[metabolite]["reactions"]:
                reactionid = self.reactions[rid]["reactionid"]
                coeffs.append(float(self.metabolites[metabolite]["reactions"][rid]))
                rids.append(self.x[reactionid])
            if output == "debug":
                print(metabolite, coeffs, rids)
            if metabolite in self.inconstant_onlyincrease:
                self.problem += pulp.lpDot(coeffs, rids) >= 0
            elif metabolite in self.inconstant:
                #self.problem += pulp.lpDot(coeffs, rids) >= 0
                pass
            else:
                self.problem += pulp.lpDot(coeffs, rids) == 0
        self.status = self.problem.solve(pulp.GLPK(msg=0))
        #self.status = self.problem.solve()

        #
        # Minimum utiization of enzymes
        #
        self.problem2 = pulp.LpProblem('sample', pulp.LpMinimize)
        self.x2 = []
        x2 = []
        #
        #for rid in self.reactionid_array:
        for i in range(len(self.reactionid_array)):
            rid = str(self.x[i])
            lb= self.reactions[rid]["lb"]
            ub= self.reactions[rid]["ub"]
            if self.x[i].value() == None:
                pass
            elif ((self.x[i].value() > 0) & (float (lb) < 0)):   #190528修正
                lb = 0.0
            elif ((self.x[i].value() < 0) & (float (ub) > 0)):   #190528修正
                ub = 0.0
            else:
                pass
            if output == "debug":
                print(rid, lb, ub)
            self.x2.append(pulp.LpVariable(rid, float(lb), float(ub), pulp.LpContinuous))
        #
        coeffs = []
        rids = []
        for i in range(len(self.reactionid_array)):
            if self.x[i].value() == None:
                coef = 0.0
            elif self.x[i].value() > 0:
                coef = 1.0
            elif self.x[i].value() < 0:
                coef = -1.0
            else:
                coef = 0.0
            coeffs.append(coef)
            rids.append(self.x2[i])
            if output == "debug":
                print(self.x2[i], coeffs, rids)
        self.problem2 += pulp.lpDot(coeffs, rids)
        #
        for rid in self.constraint:
            reactionid = self.reactions[rid]["reactionid"]
            value = self.constraint[rid]
            if output == "debug":
                print(rid, reactionid, value)
            self.problem2 +=  self.x2[reactionid] == float(value)
        for i in range(len(self.reactionid_array)):
            rid = str(self.x[i])
            if rid in self.optstep:
                value = self.x[i].value() * 0.9999
                #print(i, rid, value)
                self.problem2 +=  self.x2[i] == float(value)
        #
        for metabolite in self.maboliteid_array:
            coeffs = []
            rids = []
            for rid in self.metabolites[metabolite]["reactions"]:
                reactionid = self.reactions[rid]["reactionid"]
                coeffs.append(float(self.metabolites[metabolite]["reactions"][rid]))
                rids.append(self.x2[reactionid])
            if output == "debug":
                print(metabolite, coeffs, rids)
            if metabolite in self.inconstant_onlyincrease:
                self.problem2 += pulp.lpDot(coeffs, rids) >= 0
            elif metabolite in self.inconstant:
                #self.problem += pulp.lpDot(coeffs, rids) >= 0
                pass
            else:
                self.problem2 += pulp.lpDot(coeffs, rids) == 0
        self.status2 = self.problem2.solve(pulp.GLPK(msg=0))
        #self.status = self.problem2.solve()



        return(pulp.LpStatus[self.status], pulp.value(self.problem.objective))


    def show_reaction(self, reaction):
        """Method to show or output the optimized metabolic flux level of given reaction "reaction".
        This method can be performed after the execution of solve() or solve2().


        Args:
            reaction (str): Reaction ID

        Returns:
            nothing:

        Examples:
            >>> model.show_reaction("output.txt")
        See


        History:

        """
        print(reaction,"\t" ,self.get_value(reaction))



    def show_result2(self, filename = "none", mode = "nonzero"):
        self.show_result(filename, mode)


    def show_result(self, filename = "none", mode = "nonzero"):
        """Method to show or output the optimized metabolic flux level.
        This method can be performed after the execution of solve() or solve2().


        Args:
            filename (str): Output file is created when filename is not "none".
            mode (str): When mode is "nonzero", reactions whose flux value is within -0.0001 < value <0.0001 are not shown

        Returns:
            nothing:

        Examples:
            >>> model.show_result()
            >>> model.show_result("output.txt")
        See


        History:

        """

        if filename == "none":
            for i in range(len(self.reactionid_array)):
                rid = str(self.x2[i])
                value =  self.x2[i].value()
                group = self.reactions[rid]["group"]
                if value == None:
                    continue
                if mode == "nonzero":
                    if -0.0001 < value <0.0001:
                        continue
                print(rid,"\t" ,'{:>9.4f}'.format(value),"\t", self.reactions[rid]["reaction"],"\t", str(group))
        else:
            with open(filename, mode='w') as f:
                for i in range(len(self.reactionid_array)):
                    rid = str(self.x2[i])
                    reaction = self.reactions[rid]["reaction"]
                    lb = self.reactions[rid]["lb"]
                    ub = self.reactions[rid]["ub"]
                    group = self.reactions[rid]["group"]
                    value =  self.x2[i].value()
                    if value == None:
                        continue
                    if mode == "nonzero":
                        if -0.0001 < value <0.0001:
                            continue
                    line = rid+"\t" +str(group)+"\t" +str(lb)+"\t" +str(ub)+"\t" +str(value)+"\t"+"\t"+self.reactions[rid]["reaction"]+"\n"
                    f.writelines(line)

    def add_constraint(self, reaction, value = 0):
        """Method to constraint the metabolic flux level of "reaction" at a fixed level "value"


        Args:
            reaction (str): Reaction ID
            value (float): fixed metabolic flux level
        Returns:
            nothing

        Examples:
            >> model.add_constraint("Rxn1", value = 10.0)


        History:
        """

        self.constraint[reaction] = value

    def delete_constraint(self, reaction):
        """Method to remove a constraint of flux level from "reaction".


        Args:
            reaction (str): Reaction ID
        Returns:
            nothing

        Examples:
            >> model.delete_constraint("Rxn1")


        History:
        """

        del self.constraint[reaction]



    def get_metabolites(self):
        """Method to get a list of metaboilte ids in the metabolic model.



        Args:
            nothing:

        Returns:
            metabolitelist (list of str): list of metabolite ids in the metabolic model

        Examples:
            >>> metabolitelist = model.get_metabolites()
        See


        History:

        """

        return(self.metabolites.keys())

    def get_reactions(self):
        """Method to get a list of reaction ids in the metabolic model.



        Args:
            nothing:

        Returns:
            reactionlist (list of str): list of reaction ids in the metabolic model

        Examples:
            >>> reactionlist = model.get_reactions()
        See


        History:

        """

        return(self.reactions.keys())


    def get_reactions_for_genedeletion(self):
        """Method to obtain a list of reaction IDs (genes) available for the gene deletion simulation.

        Args:
            nothing
        Returns:
            reactions (list of str): A list of reaction IDs that are not used for the optimizaton function and constrained to a fixed value.

        Examples:
            >> reactions =  model.deletiongenes()


        History:
        """

        reactions = []
        for reaction in self.reactionid_array:
            if reaction in self.constraint:
                continue
            if reaction in self.optstep:
                continue
            reactions.append(reaction)
        return(reactions)



    def get_value(self, reaction):
        """Method to get optimized metabolic flux level of given reaction "reaction".
        This method can be performed after the execution of solve() or solve2().


        Args:
            reaction (str): Reaction ID

        Returns:
            value (float): optimized metabolic flux level of "reaction"

        Examples:
            >>> value = model.get_value("Rxn1")
        See


        History:

        """

        reactionid = self.reactions[reaction]["reactionid"]

        return(self.x2[reactionid].value())





    def wrireLog(self, filename, contentslist, output = "display"):
        """Method to output a process log to a logfile.


        Args:
            filename (str): Filename for output
            contentslist (list): the content in this list is added to the end of log file as tab-separated line
            output (str) =  Log data is no shown in display only when output = "nodisplay"

            reaction (str): Reaction ID

        Returns:
            nothing

        Examples:
            >>>model.wrireLog("logfile.txt", ["Rxn", reactionide, fluxlevel, status], output = "nodisplay" )
        See


        History:

        """

        f = open(filename, 'a',  encoding='UTF-8')
        strlist = [str(x) for x in contentslist]
        strtemp = "\t".join(strlist)
        if output != "nodisplay":
            print(strtemp)
        f.write(strtemp+"\n")
        f.close()


    def add_optstep(self, reaction, value = 1.0):
        """Method to add one reaction to be used as a optimization function.
        More than two reactions can be used as a optimization function.
        For example,

            >> model.add_optstep("rxn1", value = 1.0)
            >> model.add_optstep("rxn2, value = -0.5)

            The following objective function  F is maximized
            F = 1.0 * "rxn1" - 0.5 * "rxn2"

        Args:
            reaction (str): Reaction ID
            value (float): the optimization function is maxizied when "value" is positive value, vice versa


        Returns:
            nothing

        Examples:
            >> objfunc = "R_biomass"
            >> model.add_optstep(objfunc, value = 1.0)

        History:

        """

        self.optstep[reaction] = value
    def delete_optstep(self, reaction):

        """Method to remove one reaction from the optimization function.

        Args:
            reaction (str): Reaction ID

        Returns:
            nothing

        Examples:
            >> objfunc = "R_biomass"
            >> model.remove_optstep(objfunc, value = 1.0)

        History:

        """

        del self.optstep[reaction]
    def clear_optstep(self):
        """Method to clear optimization function.


        Args:
            reaction (str): Reaction ID


        Returns:
            nothing

        Examples:
            >> model.clear_optstep()

        History:

        """

        self.optstep = {}
    def set_boundary(self, reaction, lb, ub):
        """Method to set lower and upper boundary levels of metabolic flux level of the "reaction"


        Args:
            reaction (str): Reaction ID
            lb (float): lower boundary of metabolic flux
            ub (float): upper boundary of metabolic flux
        Returns:
            nothing

        Examples:
            >> model.set_boundary("Rxn1", 0.0, 12.0)


        History:
        """

        self.reactions[reaction]["lb"] = lb
        self.reactions[reaction]["ub"] = ub


    def set_upperboundary(self, reaction, value = 0):
        """Method to set upper boundary levels of metabolic flux level of the "reaction"


        Args:
            reaction (str): Reaction ID
            ub (float): upper boundary of metabolic flux
        Returns:
            nothing

        Examples:
            >> model.set_upperboundary("Rxn1", 12.0)


        History:
        """
        self.reactions[reaction]["ub"] = value




    def set_lowerboundary(self, reaction, value = 0):
        """Method to set lower boundary levels of metabolic flux level of the "reaction"


        Args:
            reaction (str): Reaction ID
            lb (float): lower boundary of metabolic flux
        Returns:
            nothing

        Examples:
            >> model.set_lowerboundary("Rxn1", 12.0)


        History:
        """
        self.reactions[reaction]["lb"] = value





    def add_inconstant_onlyincrease(self, metabolite):

        """Method to set a metabolite whose level is allowed to increase but not to decrease.
        Metabolites secreted from cell to medium should be set by the function.
        Ex: Lacate for human cell culture, Ethanol for yeast



        Args:
            metabolite (str): metabolite ID
        Returns:
            nothing

        Examples:
            >> model.add_inconstant_onlyincrease("glc[e]")


        History:
        """


        if type(metabolite) is list:
            for metaboliteid in metabolite:
                self.inconstant_onlyincrease[metaboliteid] = 1
        else:
            self.inconstant_onlyincrease[metabolite] = 1

    def delete_inconstant_onlyincrease(self, metabolite):
        """Method to remove a metabolite from the list of "inconstant_onlyincrease" metabolites.



        Args:
            metabolite (str): metabolite ID
        Returns:
            nothing

        Examples:
            >> delete_inconstant("glc[e]")


        History:
        """

        if type(metabolite) is list:
            for metaboliteid in metabolite:
                del self.inconstant_onlyincrease[metabolite]
        else:
            del self.inconstant_onlyincrease[metabolite]

    def add_inconstant(self, metabolite):

        """Method to set a metabolite whose level is allowed to increase and decrease.
        Metabolites incorparated into cell to medium should be set by the function.
        Ex: glucose, nitrogen source


        Args:
            metabolite (str): metabolite ID
        Returns:
            nothing

        Examples:
            >> model.add_inconstant("glc[e]")


        History:
        """



        if type(metabolite) is list:
            for metaboliteid in metabolite:
                self.inconstant[metaboliteid] = 1
        else:
            self.inconstant[metabolite] = 1
    def delete_inconstant(self, metabolite):
        """Method to remove a metabolite from the list of "inconstant" metabolites.



        Args:
            metabolite (str): metabolite ID
        Returns:
            nothing

        Examples:
            >> delete_inconstant("glc[e]")


        History:
        """

        if type(metabolite) is list:
            for metaboliteid in metabolite:
                del self.inconstant[metaboliteid]
        else:
            del self.inconstant[metabolite]


    def single_gene_deletions(self, target, filename = "output.txt"):
        """Method to perform a single gene deletion expeeriment for all available reactions


        Args:
            target (str): Reaction IDs to increase by the gene deletion
            filename (str): Name of output file name

        Returns:
            nothing

        Examples:
            >> model.single_gene_deletions("Rxn1", "sgout.txt")


        History:
        """

        reactions =  self.get_reactions_for_genedeletion()
        for reaction in  reactions:
            if reaction in self.constraint:
                continue
            self.add_constraint(reaction, value = 0)
            status = self.solve2()
            value = self.get_value(target)
            content = [reaction,value,status[0],status[1]]
            self.wrireLog(filename, content)
            self.delete_constraint(reaction)
