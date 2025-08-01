from __future__ import print_function, division
import streamlit as st
import numpy as np
import pandas as pd
import time
from hydralit import HydraHeadApp
import hydralit_components as hc
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
from itertools import chain
from Bio.pairwise2 import format_alignment
from operator import itemgetter
import os
import sys
from openpyxl import reader,load_workbook,Workbook
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode
import matplotlib.pyplot as plt
import io
from streamlit_lottie import st_lottie
import json
import requests
from pathlib import Path
from io import StringIO
import seaborn as sns
try:
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
except ImportError as exception:
    print("[!] Could not import Biopython modules", file=sys.stderr)
    raise exception
import joblib


# -------------------------------------------------------------------------------------------------------------------------------------
class PredictApp(HydraHeadApp):

    def __init__(self, title = 'Predict your peptide', **kwargs):
        # self.__dict__.update(kwargs)
        self.title = title

    def run(self):
        with open('style2.css') as f:
            st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
            with st.container():
                coll01, coll02, coll04 = st.columns((0.1,10,0.5))
                coss01, coss02, coss03 = st.columns((0.99,19,0.9))
                html_temp = """
                                                                <div style="background-color:{};height:{};width:{};">
                                                                </div>
                                                                <div style="background-color:#1F3D7C;color:white;padding:8px;border-radius:2px">
                                                                <div id="head" style="background-color:{};padding:1px;border-radius:'1px';">
                                                                </div>
                                                                """
                Ideal_datafa = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 30px; background-color: #1F3D7C; border-radius: 2px; text-align:center;"> 🧬 Predict your peptides  󠀠🧬</p>'
                coss02.markdown(Ideal_datafa, unsafe_allow_html=True)
                
                st.write("##")
                cd1,cd2, cd3 = st.columns((0.5,5,10))
                help_input="- Please enter your peptide sequence for prediction. \n- Please choose an input between FASTA text format or FASTA file format. \n- Sequence longer than 200 peptides should not be used as this may cause a delay in prediction."
                Ideala = '<div align="left"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #1F3D7C; border-radius: 5px; text-align:center;">Please enter your peptide or upload file 👇</p>'
                cd2.markdown(Ideala, unsafe_allow_html=True)
                cd3.markdown(help_input, unsafe_allow_html=True)
                clol01,clol01, clol02, clol03, clol04 = st.columns((0.5,1.5,9,0.1,0.4))
                with clol01:
                    
                    def load_lottiefile(filepath: str):
                        with open (filepath,"r") as f:
                            return json.load(f)
                    def load_lottieurl(url: str):
                        r = requests.get(url)
                        if r.status_code != 200:
                            return None
                        return r.json()
                    lottie2_codingsd = load_lottieurl("https://assets8.lottiefiles.com/packages/lf20_urdso8u9.json")
                    st_lottie(lottie2_codingsd, height=140,  key="codvingq")
                # cl1, cl2, cl3, cl4 = st.columns((0.60,15,0.009,0.5))
                
                # clol02.write("1.Enter your Peptide 👇")                                                                                                                                                                                                                                                                                                                                                    
                text_seq = clol02.text_area("",help=help_input)
                clss1, clss2, clss3, clss4,clss5 = st.columns((2,4,2.5,2,0.55))
                data_file_uploader = clss2.file_uploader('', type=['FASTA','txt'], accept_multiple_files=False)    
                if data_file_uploader is not None:
                    if data_file_uploader:
                        if (data_file_uploader.name[-5:] == ('fasta')) or (data_file_uploader.name[-5:] == ('FASTA')):
                            data_file_uploader.seek(0)
                            data_file_uploader = StringIO(data_file_uploader.getvalue().decode("utf-8"))
                            data_file = data_file_uploader.getvalue()                            
                                
                        elif (data_file_uploader.name[-3:] == 'txt') :
                            data_file_uploader = StringIO(data_file_uploader.getvalue().decode("utf-8"))
                                # st.write("filename:", data_file.name)
                            data_file = data_file_uploader.getvalue()                                   
                    else:
                        st.write('your file is incorrect')   
                
                #selection threshold 
                with clss3:
                    # st.write('Plase select range of threshold model anti/non-microbial peptide👇')
                    # st.checkbox("Disable selectbox widget", key="disabled")
                    st.markdown(
                        """
                        <style>
                        [data-baseweb="select"] {
                            margin-top: -50px;
                            background-color: #1F3D7C;
                            color: #353131;
                            border-radius: 5px;
                            font-size: 20px;
                        }
                        </style>
                        """,
                        unsafe_allow_html=True,
                    )
                    I_option_anti = '<div align="left"><p style="font-sans-serif:; color: #353131; font-size: 16px; background-color: white; border-radius: 5px; text-align:center;">Select the threshold for antimicrobial peptide model 👇</p>'
                    st.markdown(I_option_anti, unsafe_allow_html=True)
                    # st.markdown("Select the range threshold antimicrobial peptide model 👇")
                    option_anti = st.selectbox(
                        '',["50", "60", "70", "80"])


                    # st.write('You selected:', option_anti)
                    # st.write('<style>div.row-widget.stRadio > div{flex-direction:row;}</style>', unsafe_allow_html=True)
                    
                with clss4:
                    # st.checkbox("Disable selectbox widget", key="disabled")
                    st.markdown(
                        """
                        <style>
                        [data-baseweb="select"] {
                            margin-top: -50px;
                            background-color: #1F3D7C;
                            
                            color: #353131;
                            border-radius: 5px;
                            font-size: 20px;
                        }
                        </style>
                        """,
                        unsafe_allow_html=True,
                    )
                    # st.markdown("Select the range threshold gram bacteria model 👇")
                    I_option_gram = '<div align="left"><p style="font-sans-serif:; color: #353131; font-size: 16px; background-color: white; border-radius: 5px; text-align:center;">Select the threshold for <br>gram bacteria model 👇</p>'
                    st.markdown(I_option_gram, unsafe_allow_html=True)
                    option_gram = st.selectbox(
                        ' ',["50", "60", "70", "80"])

                    # st.write('You selected:', option_gram)
                    # st.write('<style>div.row-widget.stRadio > div{flex-direction:row;}</style>', unsafe_allow_html=True)
               
                cl001, ddd,cl002,cl003, cl004 = st.columns((0.95,3.8,2.5,4,16))  
                if cl002.button("🦠 Predict"):
                    try:
                        if (text_seq != '') and (data_file_uploader is not None):                                               
                            Sequence = 'two input'
                        elif (text_seq != '') and (data_file_uploader is None):
                            Sequence = text_seq
                        elif (text_seq == '') and (data_file_uploader is not None):
                            Sequence = data_file 
                        elif (text_seq == '') and (data_file_uploader is None):
                            Sequence = ''
                        if Sequence == "" :
                            cffs1,cffs2,cffs3 = st.columns((0.96,12,0.65))
                            with cffs2:
                                Ideal_er = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">**Peptide not found‼️ \n please input your peptide</p>'
                                st.markdown(Ideal_er, unsafe_allow_html=True)
                                Ideal_forexa = '<div align="left"><p style="font-sans-serif:; color: black; font-size: 15px; background-color: white; border-radius: 5px;">For example, input your peptide FASTA format</p>'
                                st.markdown(Ideal_forexa, unsafe_allow_html=True)                                  
                                html_temp = """
                                                    <div style="background-color:#D1F0FF;padding:1px">
                                                    <h8 style="color:black;text-align:left;font-size:80%;"><u>Sample1</u><br>>Sequence_name1<br>DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW<br> 
                                                    <br><u>Sample2</u><br>>Sequence_name1<br>DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW<br>>Sequence_name2<br>FPFLLSLIPSAISALKKL </h1>
                                                    </div><br>"""
                                st.markdown(html_temp,unsafe_allow_html=True)
                        elif Sequence == 'two input' :
                            st.error('Please select only one input of peptides between FASTA text format or FASTA file format.') 
                        else:
                            # model---------------------------------------------------------------------------

                            model_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "model_gbc_resize_test.joblib"))
                            # st.write("🔍 Loading model from:", model_path)
                            model_anti_or_non = joblib.load(model_path)
                            # st.success("✅ model_anti_or_non loaded")

                            model_path2 = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "model_RF_fulldata_gram-.pkl"))
                            # st.write("🔍 Loading model from:", model_path2)
                            model_angram_negative = joblib.load(model_path2)
                            # st.success("✅ model_angram_negative loaded")

                            model_path3 = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "model_RF_gramPos_resize_test.joblib"))
                            # st.write("🔍 Loading model from:", model_path3)
                            model_angram_post = joblib.load(model_path3)
                            # st.success("✅ model_angram_post loaded")

                            # input data list ----------------------------------------------------------------  
                            
                            split_sequence = Sequence.split("\n")
                            
                            list_clean_text = []
                            for i in range(len(split_sequence)):
                                clean_text = " ".join([word for word in split_sequence[i].split()])
                                list_clean_text.append(clean_text)
                                # st.write(list_clean_text)
                                    
                            seq = []
                            name_seq = []
                            y = 0
                            i = -1
                            for a in range(len(list_clean_text)):
                                if list_clean_text[a] != '':
                                    if '>' in list_clean_text[a]:
                                        name_seq.append(list_clean_text[a])
                                        i = i+1
                                    else:
                                        if len(seq) == (i):
                                            seq.append(list_clean_text[a])
                                        else:
                                            seq[i] = (str(seq[i]) + str(list_clean_text[a]))
                                        seq[i] = seq[i].replace(" ", "")
                                        seq[i] = seq[i].upper()
                            list1 = []
                            list2 = []
                                                            
                            # not_amino = ['B','J', 'O', 'U', 'X','Z']
                            list_seq =[]
                            for i in range(len(seq)):
                                if (('B') in seq[i]):
                                    Ideal_erb = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">Please check, your peptide should not have B, J, O, U, X, Z </p>'
                                    st.markdown(Ideal_erb, unsafe_allow_html=True)
                                    # st.write(seq[i]) 
                                elif (('J') in seq[i]):
                                    Ideal_erj = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">Please check, your peptide should not have B, J, O, U, X, Z </p>'
                                    st.markdown(Ideal_erj, unsafe_allow_html=True)
                                    # st.write(seq[i])
                                elif (('O') in seq[i]):
                                    Ideal_ero = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">Please check, your peptide should not have B, J, O, U, X, Z </p>'
                                    st.markdown(Ideal_ero, unsafe_allow_html=True)
                                    # st.write(seq[i])
                                elif (('U') in seq[i]):
                                    Ideal_eru = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">Please check, your peptide should not have B, J, O, U, X, Z </p>'
                                    st.markdown(Ideal_eru, unsafe_allow_html=True)
                                    # st.write(seq[i])
                                elif (('X') in seq[i]):
                                    Ideal_erx = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">Please check, your peptide should not have B, J, O, U, X, Z </p>'
                                    st.markdown(Ideal_erx, unsafe_allow_html=True)
                                    # st.write(seq[i])
                                elif (('Z') in seq[i]):
                                    Ideal_erz = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">Please check, your peptide should not have B, J, O, U, X, Z </p>'
                                    st.markdown(Ideal_erz, unsafe_allow_html=True)
                                    # st.write(seq[i])
                                else:
                                    list1.append(name_seq[i])
                                    list2.append(seq[i])                              

                            df_user_name_seq = pd.DataFrame(list(zip(list1,list2)),columns =['Name','Sequence'])
                                    
                            def CalRasidal(smi):
                                A = smi.count('A')
                                V = smi.count('V')
                                I = smi.count('I')
                                L = smi.count('L')
                                M = smi.count('M')
                                F = smi.count('F')
                                W = smi.count('W')
                                C = smi.count('C')
                                wee = A+V+I+L+M+F+C+W
                                perwee = (wee/len(smi))*100
                                hydrophobic = format(perwee, '.2f')
                                                                                
                                #hydrophilic
                                R = smi.count('R')
                                N = smi.count('N')
                                D = smi.count('D')
                                Q = smi.count('Q')
                                E = smi.count('E')
                                K = smi.count('K')
                                phii = R+N+D+Q+E+K

                                # hydrophilic.append(phii)
                                perPhi = (phii/len(smi))*100
                                hydrophilic = format(perPhi, '.2f')
                                            
                                #uncharged
                                S = smi.count('S')
                                T = smi.count('T')
                                N = smi.count('N')
                                Q = smi.count('Q')
                                too = S+T+N+Q

                                # uncharged.append(too)
                                pertoo = (too/len(smi))*100
                                uncharged = format(pertoo, '.2f')
                                            
                                #Charged
                                #% positive charge
                                K = smi.count('K')
                                R = smi.count('R')
                                H = smi.count('H')
                                Lo = K+R+H
                                perLo = (Lo/len(smi))*100
                                positiveC = format(perLo, '.2f')
                                            
                                #% Negative charge
                                D = smi.count('D')
                                E = smi.count('E')
                                so = D+E
                                NegativeC = (so/len(smi))*100
                                NegativeC = format(NegativeC, '.2f')
                                                                        
                                #Molecular Weight Calculation
                                analysed_seq = ProteinAnalysis(smi)
                                MW=  format(analysed_seq.molecular_weight(),'.2f')
                                return hydrophobic, hydrophilic, uncharged, positiveC, NegativeC, MW

                            def CalpI(smi):
                                glob = GlobalDescriptor(smi)
                                glob.isoelectric_point()
                                Po = glob.descriptor
                                Pooo = Po.tolist()
                                pI_ = list(chain.from_iterable(Pooo))
                                myList = list(np.around(np.array(pI_),2)) #ให้ pI เหลือ 2 ตำแหน่ง
                                myList = myList.pop(0)
                                strmyList = str(myList)                                  
                                return strmyList
                                            
                            def Phipho(smi):
                                num_R = smi.count('R')
                                phiValue_R = num_R*(3.0)
                                Value_R = num_R*(-2.53)
                                num_N = smi.count('N')
                                phiValue_N = num_N*(0.2)
                                Value_N = num_N*(-0.78)

                                num_D = smi.count('D')
                                phiValue_D = num_D*(3.0)
                                Value_D = num_D*(-0.90)

                                num_Q = smi.count('Q')
                                phiValue_Q = num_Q*(0.2)
                                Value_Q = num_Q*(-0.85)

                                num_E = smi.count('E')
                                phiValue_E = num_E*(3.0)
                                Value_E = num_E*(-0.74)

                                num_K = smi.count('K')
                                phiValue_K = num_K*(3.0)
                                Value_K = num_K*(-1.5)
                                num_H = smi.count('H')
                                phiValue_H = num_H*(-0.5)
                                Value_H = num_H*(-0.40)

                                num_S = smi.count('S')
                                phiValue_S = num_S*(0.3)
                                Value_S = num_S*(-0.18)

                                num_T = smi.count('T')
                                phiValue_T = num_T*(-0.4)
                                Value_T = num_T*(-0.05)

                                num_A = smi.count('A')
                                phiValue_A = num_A*(-0.5)
                                Value_A = num_A*(0.62)

                                num_C = smi.count('C')
                                phiValue_C = num_C*(-0.1)
                                Value_C = num_C*(0.29)
                                num_F = smi.count('F')
                                phiValue_F = num_F*(-2.5)
                                Value_F = num_F*(1.19)

                                num_G = smi.count('G')
                                phiValue_G = num_G*(0)
                                Value_G = num_G*(0.48)

                                num_I = smi.count('I')
                                phiValue_I = num_I*(-1.8)
                                Value_I = num_I*(1.38)

                                num_L = smi.count('L')
                                phiValue_L = num_L*(-1.8)
                                Value_L = num_L*(1.06)

                                num_M = smi.count('M')
                                phiValue_M = num_M*(-1.3)
                                Value_M = num_M*(0.64)

                                num_P = smi.count('P')
                                phiValue_P = num_P*(0)
                                Value_P = num_P*(0.12)
                                num_V = smi.count('V')
                                phiValue_V = num_V*(-1.5)
                                Value_V = num_V*(1.08)

                                num_W = smi.count('W')
                                phiValue_W = num_W*(-3.4)
                                Value_W = num_W*(0.81)

                                num_Y = smi.count('Y')
                                phiValue_Y = num_Y*(-2.3)
                                Value_Y = num_Y*(0.26)

                                #listPhi_value = [Value_D, Value_E, Value_K, Value_N, Value_Q, Value_R]
                                list_phivalue = [phiValue_D, phiValue_E, phiValue_K, phiValue_N, phiValue_Q, phiValue_R, phiValue_A, phiValue_C, phiValue_F, phiValue_G, phiValue_I, phiValue_L, phiValue_M, phiValue_P, phiValue_V, phiValue_W, phiValue_Y]
                                list_value = [Value_D, Value_E, Value_K, Value_N, Value_Q, Value_R,Value_A, Value_C, Value_F, Value_G, Value_I, Value_L, Value_M, Value_P, Value_V, Value_W, Value_Y]

                                score_hydrophilic = format(sum(list_phivalue)/len(smi),'.2f')
                                Score_hydrophobic = format(sum(list_value)/len(smi),'.2f')
                                # st.write('gggg'+ score_hydrophilic[0])
                                # print('score hydrophilic:',Sum1,'Score hydrophobic:',Sum2)
                                return score_hydrophilic,Score_hydrophobic
                                        
                                            
                            Hydrophobic_list =['F','W','I','L','V','A','M','C']
                            Neutral_list = ['Y','H','T','S','P','G']
                            Hydrophilic_list =['R','Q','K','N','E','D']
                            charged_list = ['K','R','H','E','D']
                            Uncharged_list = ['T','S','Q','N']
                            sim_list = ['DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW', 'MKFTIVFLLLACVFAMGVATPGKPRPYSPRPTSHPRPIRVRREALAIEDHLTQAAIRPPPILPA', 'MASTERNFLLLSLVVSALSGLVHRSDAAEISFGSCTPQQSDERGQCVHITSCPYLANLLMVEPKTPAQRILLSKSQCGLDNRVEGLVNRILVCCPQSMRGNIMDSEPTPSTRDALQQGDVLPGNDVCGFLFADRIFGGTNTTLWEFPWMVLLQYKKLFSETYTFNCGGALLNSRYVLTAGHCLASRELDKSGAVLHSVRLGEWDTRTDPDCTTQMNGQRICAPKHIDIEVEKGIIHEMYAPNSVDQRNDIALVRLKRIVSYTDYVRPICLPTDGLVQNNFVDYGMDVAGWGLTENMQPSAIKLKITVNVWNLTSCQEKYSSFKVKLDDSQMCAGGQLGVDTCGGDSGGPLMVPISTGGRDVFYIAGVTSYGTKPCGLKGWPGVYTRTGAFIDWIKQKLEP','EDLTVKIGDFGLATEKSRWSGSHQFEQLS','MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR', 'MTTCSRQFTSSSSMKGSCGIGGGIGAGSSRISSVLAGGSCRAPNTYGGGLSVSSSRFSSGGAYGLGGGYGGGFSSSSSSFGSGFGGGYGGGLGAGLGGGFGGGFAGGDGLLVGSEKVTMQNLNDRLASYLDKVRALEEANADLEVKIRDWYQRQRPAEIKDYSPYFKTIEDLRNKILTATVDNANVLLQIDNARLAADDFRTKYETELNLRMSVEADINGLRRVLDELTLARADLEMQIESLKEELAYLKKNHEEEMNALRGQVGGDVNVEMDAAPGVDLSRILNEMRDQYEKMAEKNRKDAEEWFFTKTEELNREVATNSELVQSGKSEISELRRTMQNLEIELQSQLSMKASLENSLEETKGRYCMQLAQIQEMIGSVEEQLAQLRCEMEQQNQEYKILLDVKTRLEQEIATYRRLLEGEDAHLSSSQFSSGSQSSRDVTSSSRQIRTKVMDVHDGKVVSTHEQVLRTKN']

                            def calculate_align_sequences(sequence_A, sequence_B, **kwargs):
                                def _calculate_identity(sequenceA, sequenceB):
                                    sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
                                    c_matches = [sa[i] == sb[i] for i in range(sl)]
                                    matches = sum([1 for i in range(sl) if (sa[i] == sb[i])])
                                    # gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
                                    # gap_id = (100 * sum(matches)) / gapless_sl
                                    gaps_al = sum([1 for i in range(sl) if (sa[i] == '-' or sb[i] == '-')])
                                    return (matches, gaps_al)
                                            
                                def _calculate_similarity(sequenceA, sequenceB):
                                    sa, sb, sal, sbl = sequenceA, sequenceB, len(sequenceA), len(sequenceB)
                                    c_sim = 0
                                    align = []
                                    for i in range(sal):
                                        if ((sa[i] in Hydrophobic_list and sb[i] in Hydrophobic_list ) or (sa[i] in Neutral_list and sb[i] in Neutral_list ) or (sa[i] in Hydrophilic_list and sb[i] in Hydrophilic_list) or (sa[i] in charged_list and sb[i] in charged_list)  or (sa[i] in Uncharged_list and sb[i] in Uncharged_list) or (sa[i] == sb[i])):
                                        # if ((sa[i] in Aliphatic_list and sb[i] in Aliphatic_list ) or (sa[i] in Hydrophobic_list and sb[i] in Hydrophobic_list ) or (sa[i] in Hydroxyl and sb[i] in Hydroxyl) or (sa[i] in charged_list and sb[i] in charged_list) or (sa[i] in Acidic and sb[i] in Acidic) or (sa[i] in Hydrophilic_list and sb[i] in Hydrophilic_list) or (sa[i] in a_list and sb[i] in a_list) or (sa[i] in b_list and sb[i] in b_list) or (sa[i] in c_list and sb[i] in c_list) or (sa[i] == sb[i]) or (sa[i] in a and sb[i] in a)):
                                            if (sa[i] != sb[i]):
                                                align.append(':')
                                            elif (sa[i] == sb[i]):
                                                align.append('|')
                                            elif ((sa[i] != sb[i]) and (sa[i] != '-') and (sb[i] != '-')):
                                                align.append('.')
                                            c_sim = c_sim+1
                                        else:
                                            if ((sa[i] != sb[i]) and (sa[i] != '-') and (sb[i] != '-')):
                                                align.append('.')
                                            else:  
                                                align.append(' ')
                                    align_l = ''.join(align)
                                    return (align_l, c_sim)
                                        
                                matrix = kwargs.get('matrix', matlist.blosum62)
                                gap_open = kwargs.get('gap_open', -10.0)
                                gap_extend = kwargs.get('gap_extend', -0.5)

                                alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                                                matrix, gap_open, gap_extend,
                                                                penalize_end_gaps=(False, False) )
                                seq_id_list = []
                                x = 0
                                        
                                for alignment in alns: 
                                    aligned_A, aligned_B, score, begin, end = alignment
                                    seq_id = _calculate_identity(aligned_A, aligned_B)
                                    seq_id_list.append(seq_id)
                                inden = seq_id_list[0]
                                for i in seq_id_list:
                                    if i <= inden:
                                        aligned_A, aligned_B, score, begin, end = alns[x]
                                    x=x+1
                                            
                                # Calculate sequence identity -----------------------------------------------------------------
                                matches, gaps_al = _calculate_identity(aligned_A, aligned_B)
                                align_l, c_sim  = _calculate_similarity(aligned_A, aligned_B)

                                c_dot = sum([1 for i in range(len(align_l)) if (align_l[i] == '.')])
                                c_co = sum([1 for i in range(len(align_l)) if (align_l[i] == ':')])
                                len_a, len_b = len(sequence_A), len(sequence_B)
                                len_al = (len_a + len_b) - (matches + c_dot + c_co)

                                gaps = (100 * gaps_al)/ len_al
                                sim = (100 * c_sim) / len_al
                                iden = (100 * matches) / len_al

                                # return (aligned_A,align_list, aligned_B), 'score: %.1f' %score, 'identity: %.1f' %seq_id,'similarity: %.1f' %sim_p
                                return (aligned_A,align_l, aligned_B), float(f'{sim:,.2f}'), float(f'{iden:,.2f}'), float(f'{gaps:,.2f}'), int(c_sim), int(matches), int(gaps_al), int(len_al)

                            def align_sequences(sequence_B):
                                # sim_list = ['DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW', 'MKFTIVFLLLACVFAMGVATPGKPRPYSPRPTSHPRPIRVRREALAIEDHLTQAAIRPPPILPA', 'MASTERNFLLLSLVVSALSGLVHRSDAAEISFGSCTPQQSDERGQCVHITSCPYLANLLMVEPKTPAQRILLSKSQCGLDNRVEGLVNRILVCCPQSMRGNIMDSEPTPSTRDALQQGDVLPGNDVCGFLFADRIFGGTNTTLWEFPWMVLLQYKKLFSETYTFNCGGALLNSRYVLTAGHCLASRELDKSGAVLHSVRLGEWDTRTDPDCTTQMNGQRICAPKHIDIEVEKGIIHEMYAPNSVDQRNDIALVRLKRIVSYTDYVRPICLPTDGLVQNNFVDYGMDVAGWGLTENMQPSAIKLKITVNVWNLTSCQEKYSSFKVKLDDSQMCAGGQLGVDTCGGDSGGPLMVPISTGGRDVFYIAGVTSYGTKPCGLKGWPGVYTRTGAFIDWIKQKLEP','EDLTVKIGDFGLATEKSRWSGSHQFEQLS','MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR', 'MTTCSRQFTSSSSMKGSCGIGGGIGAGSSRISSVLAGGSCRAPNTYGGGLSVSSSRFSSGGAYGLGGGYGGGFSSSSSSFGSGFGGGYGGGLGAGLGGGFGGGFAGGDGLLVGSEKVTMQNLNDRLASYLDKVRALEEANADLEVKIRDWYQRQRPAEIKDYSPYFKTIEDLRNKILTATVDNANVLLQIDNARLAADDFRTKYETELNLRMSVEADINGLRRVLDELTLARADLEMQIESLKEELAYLKKNHEEEMNALRGQVGGDVNVEMDAAPGVDLSRILNEMRDQYEKMAEKNRKDAEEWFFTKTEELNREVATNSELVQSGKSEISELRRTMQNLEIELQSQLSMKASLENSLEETKGRYCMQLAQIQEMIGSVEEQLAQLRCEMEQQNQEYKILLDVKTRLEQEIATYRRLLEGEDAHLSSSQFSSGSQSSRDVTSSSRQIRTKVMDVHDGKVVSTHEQVLRTKN']
                                c = 0
                                list_sim_align = []
                                list_iden_align = []
                                list_gaps_align = []
                                list_c_sim =[]
                                list_matches = []
                                list_gaps_al = []
                                list_len_al = []
                                for i in sim_list:
                                    align = calculate_align_sequences(i, sequence_B)
                                    list_sim_align.append(align[1])
                                    list_iden_align.append(align[2])
                                    list_gaps_align.append(align[3])
                                    list_c_sim.append(align[4])
                                    list_matches.append(align[5])
                                    list_gaps_al.append(align[6])
                                    list_len_al.append(align[7])
                                    
                                    #ค่าในlist --> defensin,Drosocin,Spaetzle,B-RAF,hemoglobin,keratin 1=%sim, 2=iden, 3=gaps, 4 = c_sim, 5 = matches, 6 = gaps_al, 7 = len_al
                                return list_sim_align, list_iden_align, list_gaps_align, list_c_sim, list_matches, list_gaps_al, list_len_al
                                    
                            def all_data_user(Sequence_len, hydrophobic, hydrophilic, uncharged, positive_charge, Negative_charge, Molecular_Weight, pI, score_hydrophilic, Score_hydrophobic, similarity_Betadefensin, similarity_Drosocin, similarity_Spaetzle, similarity_BRAF, similarity_hemoglobin, similarity_keratin):
                                dict = {'Sequence_len': Sequence_len, '%Hydrophobic': hydrophobic, '%Hydrophilic': hydrophilic, 
                                                '%Uncharged':uncharged, '%Positive_Charge' : positive_charge, '%Negative_Charge' : Negative_charge, 
                                                'Molecular_Weight': Molecular_Weight, 'Isoelectric_Point': pI, 'Score_Hydrophilic' :score_hydrophilic, 
                                                'Score_Hydrophobic':Score_hydrophobic, 'Similarity_Beta-defensin_1': similarity_Betadefensin,
                                                'Similarity_Drosocin': similarity_Drosocin, 'Similarity_Spaetzle': similarity_Spaetzle, 
                                                'Similarity_B-RAF' : similarity_BRAF, 'Similarity_Hemoglobin': similarity_hemoglobin, 'Similarity_Keratin' : similarity_keratin} 
                                
                                df_use_in_model = pd.DataFrame(dict, dtype = float)
                                
                                return df_use_in_model
                            # function model for predict peptide ---------------------------------------------------------------- 
                            def use_model(data_user_features_user_in_model, data_user_nec_pos_in_model):
                                # progress bar
                                progress_text = "Operation in progress. Please wait."
                                my_bar = st.progress(0, text=progress_text)
                                for percent_complete in range(100):
                                    time.sleep(0.01)
                                    my_bar.progress(percent_complete*100/len(data_user_features_user_in_model), text=progress_text)
                                    
                                # list_test_nom = data_user_features_user_in_model.values.tolist()
                                # predictions_anti_or_non = model_anti_or_non.predict(data_user_features_user_in_model)
                                real_probs_anti_or_non = model_anti_or_non.predict_proba(data_user_features_user_in_model)[0]
                                # predictions_nec = model_angram_negative.predict(data_user_nec_pos_in_model)
                                real_probs_nec = model_angram_negative.predict_proba(data_user_nec_pos_in_model)[0]
                                # predictions_pos = model_angram_post.predict(data_user_nec_pos_in_model)
                                real_probs_pos = model_angram_post.predict_proba(data_user_nec_pos_in_model)[0]
                                probs_anti_or_non = [np.round(x,5) for x in real_probs_anti_or_non]
                                probs_nec = [np.round(x,5) for x in real_probs_nec]
                                probs_pos = [np.round(x,5) for x in real_probs_pos]                                                               
                                if probs_anti_or_non[1] >= (int(option_anti)/100):                                        
                                    anti_or_non.append('antimicrobial')                               
                                    probs_anti_or_non_list.append(probs_anti_or_non[1]) 
                                    # st.subheader('✔️ Your peptide is an antimicrobial peptide.')
                                    # st.text('Probability is '+ str((probs_anti_or_non)[1]))
                                    if (probs_pos[1] >= (int(option_gram)/100)) and (real_probs_nec[1] >= (int(option_gram)/100)):
                                            # st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-positive ✚ bacteria.')
                                            # st.text('Probability is '+ str((probs_pos)[1]))
                                            # st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-negative ▬ bacteria.')
                                            # st.text('Probability is '+ str((probs_nec)[1]))
                                        pos_ro_nec.append('gram+,gram-')
                                        probs_nec_list.append(probs_pos[1]) 
                                        probs_poe_list.append(probs_nec[1])
                                    elif (real_probs_pos[1] >= (int(option_gram)/100)) and (real_probs_nec[1] < (int(option_gram)/100)):
                                        # st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-positive ✚ bacteria.')
                                        # st.text('Probability is '+ str((probs_pos)[1]))
                                        pos_ro_nec.append('gram+') 
                                        probs_poe_list.append((probs_pos)[1])
                                        probs_nec_list.append('-')
                                    elif (real_probs_pos[1] < (int(option_gram)/100)) and (real_probs_nec[1] >= (int(option_gram)/100)):
                                            # st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-negative ▬ bacteria.' )
                                            # st.text('Probability is '+ str((probs_nec)[1]))
                                        pos_ro_nec.append('gram-')
                                        probs_poe_list.append('-')   
                                        probs_nec_list.append((probs_nec)[1])   
                                    elif (real_probs_pos[1] < (int(option_gram)/100)) and (real_probs_nec[1] < (int(option_gram)/100)):
                                            # st.success(' 󠀠 󠀠✔️ 󠀠 Resist other gram of bacteria.')
                                        pos_ro_nec.append('other gram')
                                        probs_nec_list.append((probs_nec)[1]) 
                                        probs_poe_list.append((probs_pos)[1])                                                                                
                                        
                                elif real_probs_anti_or_non[1] < (int(option_anti)/100):                                        
                                    anti_or_non.append('non antimicrobial')
                                    pos_ro_nec.append("-")
                                    probs_anti_or_non_list.append(probs_anti_or_non[1])
                                    probs_nec_list.append('-') 
                                    probs_poe_list.append('-')
                                        # st.subheader('❌ Your peptide is non antimicrobial peptide.')
                                        # st.text('Probability is '+ str((probs_anti_or_non)[0]))
                                                                    
                                return anti_or_non, pos_ro_nec, probs_anti_or_non_list, probs_nec_list, probs_poe_list
        
                            len_list = []                        
                            hydrophobic_list = []
                            hydrophilic_list = [] 
                            uncharged_list = []
                            positive_charge_list = []
                            Negative_charge_list = []
                            Molecular_Weight_list = []
                            pI_list = []
                            score_hydrophilic_list = []
                            Score_hydrophobic_list = []
                            similarity_Betadefensin = []
                            similarity_Drosocin = []
                            similarity_Spaetzle = []
                            similarity_BRAF = []
                            similarity_hemoglobin = []
                            similarity_keratin = []
                            anti_or_non = []
                            pos_ro_nec = []
                            probs_anti_or_non_list = []
                            probs_nec_list = [] 
                            probs_poe_list =[]   
                            
                            # -- สี progress bar --------------------------------------------------------------
                            st.markdown("""
                                <style>
                                .stProgress > div > div > div > div {
                                    background-color: #f48fb1;
                                }
                                </style>
                            """, unsafe_allow_html=True)
                            # ---- progress bar --------------------------------------------------------------
                            progress_bar = st.progress(0)
                            status_text = st.empty()
                            total = len(df_user_name_seq)
                            # ------------------------------------------------------------------------------
                            for i in df_user_name_seq['Sequence']:

                                len_list.append(len(i))

                                hydrophobic, hydrophilic, uncharged, positiveC, NegativeC, MW = CalRasidal(i)
                                hydrophobic_list.append(hydrophobic)
                                hydrophilic_list.append(hydrophilic)
                                uncharged_list.append(uncharged)
                                positive_charge_list.append(positiveC)
                                Negative_charge_list.append(NegativeC)
                                Molecular_Weight_list.append(MW)

                                pI_pi = CalpI(i)
                                pI_list.append(pI_pi)

                                score_hydrophilic, Score_hydrophobic = Phipho(i)
                                score_hydrophilic_list.append(score_hydrophilic)
                                Score_hydrophobic_list.append(Score_hydrophobic)
                                        
                                list_sim_align, list_iden_align, list_gaps_align, list_c_sim, list_matches, list_gaps_al, list_len_al = align_sequences(i)
                                similarity_Betadefensin.append(list_sim_align[0])
                                similarity_Drosocin.append(list_sim_align[1])
                                similarity_Spaetzle.append(list_sim_align[2])
                                similarity_BRAF.append(list_sim_align[3])
                                similarity_hemoglobin.append(list_sim_align[4])
                                similarity_keratin.append(list_sim_align[5])
                                                        
                            df_use_in_model = all_data_user(len_list, hydrophobic_list, hydrophilic_list, uncharged_list, positive_charge_list, Negative_charge_list, Molecular_Weight_list, pI_list, score_hydrophilic_list, Score_hydrophobic_list, similarity_Betadefensin, similarity_Drosocin, similarity_Spaetzle, similarity_BRAF, similarity_hemoglobin, similarity_keratin)
                            
                            list_mean_ant_non = [95.98, 41.4, 33.1, 16.49, 20.49, 7.91, 10804.93, 8.79, 0.09, 0.12, 5.17, 5.03, 3.44, 5.57, 6.63, 4.41]
                            list_std_ant_non = [108.19, 11.88, 11.91, 9.65, 12.68, 7.25, 12119.58, 2.72, 0.46, 0.36, 5.08, 4.26, 3.67, 4.46, 5.05, 4.68]

                            list_mean_nor = [27.15, 44.82, 33.16, 11.2, 27.09, 3.47, 3105.05, 10.39, 0.07, 0.09, 13.08, 9.97, 3.35, 14.35, 6.78, 2.7]
                            list_std_nor = [25.87, 14.33, 14.54, 9.59, 14.83, 6.04, 2849.90, 2.02, 0.57, 0.47, 11.23, 8.3, 2.75, 10.1, 4.82, 2.11]
                            
                            df_ant_non_normed = (df_use_in_model.sub(list_mean_ant_non, axis='columns')).div(list_std_ant_non)
                            df_pos_nec_normed = (df_use_in_model.sub(list_mean_nor, axis='columns')).div(list_std_nor)                        
                            
                            # with st.container():
                            with open('style2.css') as f:
                                st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
                                        
                                # <h1 style="color:{};text-align:center;">Dataframe of predict your peptide</h1>
                                l_col1, l_col2, lasti = st.columns((0.60,12,0.6))
                                
                                with l_col2:
                                    if len(df_ant_non_normed) <= 50:
                                        for i in range(len(df_ant_non_normed)):                                   
                                        
                                            with open('style2.css') as f:
                                                with st.expander('Describe detail information'):
                                                    st.info(
                                                            """ 
                                                                - Name = Name of your sequence peptide.
                                                                - Probability = The probability that your peptide is antimicrobial peptide.
                                                                - Feature = Feature of your sequence peptide.
                                                                - Similarity = The Similarity of your peptide compare with another peptide such as antimicrobial peptide e.g. Defensin, Drosocin, Spaetzle, or non-antimicrobial e.g. B-RAF, Hemoglobin, Keratin.
                                                                - Identity = The distinguish character of your peptide compared with another peptide.
                                                                - Gaps = A break or space in the peptide compared with other peptide
                                                                """)
                                                st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
                                                
                                                cc1,cc2,cc3,cc4,cc5,cc6 = st.columns((4,4,3,3,3,3))

                                                cc1.write("""<style>.font-family: Poppins, sans-serif; {font-size:15px !important;}</style>""", unsafe_allow_html=True)
                                                cc1.write('👉🏻 Name: '+ df_user_name_seq['Name'][i])
                                                cc1.write('👉🏻 Sequence: '+ df_user_name_seq['Sequence'][i])
                                                
                                                # Predict result form model
                                                html_temp = """
                                                            <div style="background-color:{};height:{};width:{};">
                                                            </div>
                                                            <div style="background-color:#1F3D7C;color:white;padding:2px;border-radius:5px">
                                                            <div id="head" style="background-color:{};padding:2px;border-radius:'5px';">
                                                            </div>
                                                            """
                                                
                                            with cc2: 
                                                # st.write("##")
                                                # st.markdown(html_temp.format('white','2px', '60%','#1F3D7C','#FFFFFF'),unsafe_allow_html=True)

                                                anti_or_non, pos_ro_nec, probs_anti_or_non_list, probs_nec_list, probs_poe_list = use_model(df_ant_non_normed.iloc[[i]], df_pos_nec_normed.iloc[[i]])
                                                
                                                if anti_or_non[i] == 'antimicrobial':
                                                   
                                                    # st.code('✔️ 󠀠'+ anti_or_non[i])
                                                    # st.write('Probability is:'+ " 󠀠 󠀠 󠀠 "+ str(probs_anti_or_non_list[i]))
                                                    # st.code('Active against')
                                                    # Ideala = '<div align="left"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #1F3D7C; border-radius: 5px; text-align:center;">Please enter your peptide or File upload 👇</p>'
                                                    potential_anti = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:left;"> 󠀠 󠀠 Probability: 󠀠 󠀠 󠀠 𝒀𝒆𝒔✔️</p>'
                                                    st.markdown(potential_anti, unsafe_allow_html=True)
                                                    # st.markdown('Potential to be AMPs:'+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 󠀠✔️')
                                                    st.write('Probability:'+ " 󠀠 󠀠 󠀠 " + str(probs_anti_or_non_list[i]))
                                                    potential_targ = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:left;"> 󠀠 󠀠 Target Bacteria</p>'
                                                    st.markdown(potential_targ, unsafe_allow_html=True)
                                                    # st.code('Target Bacteria')

                                                    if pos_ro_nec[i] == 'gram+,gram-':
                                                        st.write('Potential againt Gram + Bacteria:'+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 󠀠✔️'+ " 󠀠 󠀠 󠀠 "+ 'Probability:'+ " 󠀠 󠀠 󠀠 " + str(probs_poe_list[i]))
                                                        # st.write('Probability:'+ " 󠀠 󠀠 󠀠 " + str(probs_poe_list[i]))
                                                        # st.write("gram+"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_poe_list[i]))
                                                        st.write('Potential againt Gram - Bacteria:'+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 ✔️'+ " 󠀠 󠀠 󠀠 "+ 'Probability:'+ " 󠀠 󠀠 󠀠 " + str(probs_nec_list[i]))
                                                        # st.write('Probability:'+ " 󠀠 󠀠 󠀠 " + str(probs_poe_list[i]))
                                                        # st.write("gram-"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_nec_list[i]))
                                                    elif pos_ro_nec[i] == 'gram+':
                                                        st.write("Potential againt Gram + Bacteria:"+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 ✔️'+ " 󠀠 󠀠 󠀠 "+'Probability:'+ " 󠀠 󠀠 󠀠 "+str(probs_poe_list[i]))
                                                    elif pos_ro_nec[i] == 'gram-':
                                                        st.write("Potential againt Gram - Bacteria:"+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 ✔️'+ " 󠀠 󠀠 󠀠 "+'Probability:'+ " 󠀠 󠀠 󠀠 "+str(probs_nec_list[i]))
                                                elif anti_or_non[i] == 'non antimicrobial':
                                                    # st.markdown('Potential to be AMPs:'+ " 󠀠 󠀠 󠀠 "+'𝑵𝒐❌ 󠀠')
                                                    potential_non = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:left;"> 󠀠 Probability:  󠀠 󠀠 󠀠 𝑵𝒐 ❌ 󠀠 </p>'
                                                    st.markdown(potential_non, unsafe_allow_html=True)
                                                    st.write('Probability:'+ " 󠀠 󠀠 󠀠 "+ str(probs_anti_or_non_list[i])) 
                                                    
                                            
                                                # if pos_ro_nec[i] == 'gram+,gram-':
                                                #     st.write("gram+"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_poe_list[i]))
                                                #     st.write("gram-"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_nec_list[i]))
                                                # elif pos_ro_nec[i] == 'gram+':
                                                #     st.write("gram+"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_poe_list[i]))
                                                # elif pos_ro_nec[i] == 'gram-':
                                                #     st.write("gram-"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_nec_list[i]))
                                                
                                                
                                                # st.write(probs_nec_list[i])
                                                # st.write(probs_poe_list[i])
                                                # st.markdown(html_temp.format('white','1px', '60%','#1F3D7C','#FFFFFF'),unsafe_allow_html=True)   
                                        
                                            # display result --------------------------------------------------------------------------
                                            # with st.container():
                                            #     left_col1, left_col2, center, right_col1, right_col2,last = st.columns((0.50,5.5,8,8,8,0.3))
                                            with cc3:
                                                # help_feature="Feature of your peptide \n"
                                                # with st.expander("Feature"):
                                                #     st.info('The identity of your peptide compare identity with another peptide.')
                                                # st.code('Feature')
                                                mark_Feature = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:center;"> 󠀠 󠀠Feature 󠀠 </p>'
                                                st.markdown(mark_Feature, unsafe_allow_html=True)
                                                #show Sequence len-------------------------------
                                                st.write('🔻 Sequence len: '+ str(len_list[i]), unsafe_allow_html=True)

                                                        #show Hydrophobic-------------------------------
                                                st.write('🔻 Hydrophobic: '+ (hydrophobic_list[i]), unsafe_allow_html=True)

                                                        #show Hydrophilic-------------------------------
                                                st.write('🔻 Hydrophilic: '+ (hydrophilic_list[i]), unsafe_allow_html=True)
                                            
                                            with cc4:
                                                # st.code('Similarity Comparison')
                                                mark_Similarity = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:center;"> 󠀠 Similarity Comparison 󠀠 </p>'
                                                st.markdown(mark_Similarity, unsafe_allow_html=True)
                                                st.write('Defensin is: '+ str(list_sim_align[0])+"%")    
                                                st.write('Drosocin is: '+ str(list_sim_align[1])+"%")
                                                st.write('Spaetzle is: '+ str(list_sim_align[2])+"%")
                                            with cc5:
                                                # st.code('Identity Comparison')
                                                mark_Identity = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:center;"> 󠀠 Identity Comparison 󠀠 </p>'
                                                st.markdown(mark_Identity, unsafe_allow_html=True)
                                                st.write('Defensin is: '+ str(list_iden_align[0])+"%")    
                                                st.write('Drosocin is: '+ str(list_iden_align[1])+"%")
                                                st.write('Spaetzle is: '+ str(list_iden_align[2])+"%")
                                            with cc6:
                                                # st.code('Gaps Comparison')
                                                mark_Gaps = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:center;"> 󠀠 Gaps Comparison 󠀠 </p>'
                                                st.markdown(mark_Gaps, unsafe_allow_html=True)
                                                st.write('Defensin is: '+ str(list_gaps_align[0])+"%")    
                                                st.write('Drosocin is: '+ str(list_gaps_align[1])+"%")
                                                st.write('Spaetzle is: '+ str(list_gaps_align[2])+"%")

                                            with open('style2.css') as f:
                                                st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
                                                with st.expander('🔵 󠀠 󠀠Select for more detail...'):
                                                    with open('style2.css') as f:
                                                        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
                                                        ff1,ff2,ff4 = st.columns((8,9,1))
                                                        with ff2:
                                                            Ideal_Feature = '<div align="center"><p style="font-family:; color:#1F3D7C; font-size: 18px; background-color: #C2DFFF;">Feature</p>'
                                                            st.markdown(Ideal_Feature, unsafe_allow_html=True)
                                                    cooll1, cooll2, cooll3,cooll4,cooll5,cooll6 = st.columns((0.60,6.5,3.5,0.5,3.5,1.1))
                                                    with cooll2:
                                                        st.write('👉🏻 Name: '+ df_user_name_seq['Name'][i])
                                                        st.write('👉🏻 Sequence: '+ df_user_name_seq['Sequence'][i])
                                                        if anti_or_non[i] == 'antimicrobial':
                                                            # st.code('✔️ 󠀠'+ anti_or_non[i])
                                                            # st.write('Probability is:'+ " 󠀠 󠀠 󠀠 "+ str(probs_anti_or_non_list[i]))
                                                            # st.code('Active against')
                                                            potential_anti2 = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:left;"> 󠀠 󠀠 Probability: 󠀠 󠀠 󠀠 𝒀𝒆𝒔✔️</p>'
                                                            st.markdown(potential_anti2, unsafe_allow_html=True)
                                                            st.write('A probability threshold of:'+ " 󠀠 󠀠 󠀠 " + option_anti + " 󠀠 󠀠 󠀠 " + " is "+ " 󠀠 󠀠 󠀠 " +str(probs_anti_or_non_list[i]))
                                                            potential_targ2 = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:left;"> 󠀠 󠀠 Target Bacteria</p>'
                                                            st.markdown(potential_targ2, unsafe_allow_html=True)
                                                            st.write('A probability threshold of:'+ " 󠀠 󠀠 󠀠 " + option_gram)

                                                            if pos_ro_nec[i] == 'gram+,gram-':
                                                                st.write('Potential againt Gram + Bacteria:'+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 󠀠✔️'+ " 󠀠 󠀠 󠀠 "+ 'Probability:'+ " 󠀠 󠀠 󠀠 " + str(probs_poe_list[i]))
                                                                st.write('Potential againt Gram - Bacteria:'+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 ✔️'+ " 󠀠 󠀠 󠀠 "+ 'Probability:'+ " 󠀠 󠀠 󠀠 " + str(probs_nec_list[i]))
                                                                # st.write("✔️ 󠀠gram+"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_poe_list[i]))
                                                                # st.write("✔️ 󠀠gram-"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_nec_list[i]))
                                                            elif pos_ro_nec[i] == 'gram+':
                                                                # st.write("✔️ 󠀠gram+"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_poe_list[i]))
                                                                st.write("Potential againt Gram + Bacteria:"+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 󠀠✔️'+'Probability:'+ " 󠀠 󠀠 󠀠 "+str(probs_poe_list[i]))
                                                            elif pos_ro_nec[i] == 'gram-':
                                                                # st.write("✔️ 󠀠gram-"+ " 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠 󠀠"+'Probability is:'+ " 󠀠 󠀠 󠀠 "+str(probs_nec_list[i]))
                                                                st.write("Potential againt Gram - Bacteria:"+ " 󠀠 󠀠 󠀠 " + ' 𝒀𝒆𝒔 󠀠✔️'+'Probability:'+ " 󠀠 󠀠 󠀠 "+str(probs_nec_list[i]))
                                                        elif anti_or_non[i] == 'non antimicrobial':
                                                            # st.code('❌ 󠀠'+ anti_or_non[i])
                                                            # st.write('Probability is:'+ " 󠀠 󠀠 󠀠 "+ str(probs_anti_or_non_list[i]))
                                                            potential_non2 = '<div align="center"><p style="font-sans-serif:; color:white; font-size: 16px; background-color: #1F3D7C; border: 2px solid #06BBCC; border-radius: 5px; text-align:left;"> 󠀠 󠀠Probability:  󠀠 󠀠 󠀠 𝑵𝒐 ❌ 󠀠 </p>'
                                                            st.markdown(potential_non2, unsafe_allow_html=True)
                                                            st.write('A probability threshold of:'+ " 󠀠 󠀠 󠀠 " + option_anti + " 󠀠 󠀠 󠀠 " + " is "+ str(probs_anti_or_non_list[i])) 
                                                        # if anti_or_non[i] == "antimicrobial":
                                                        #     st.subheader('✔️ Your peptide is an antimicrobial peptide.')
                                                        #     st.text('Probability is '+ str((probs_anti_or_non_list)[i]))
                                                        #     if (pos_ro_nec[i] == "gram+,gram-"):
                                                        #         st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-positive ✚ bacteria.')
                                                        #         st.text('Probability is '+ str((probs_poe_list)[i]))
                                                        #         st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-negative ▬ bacteria.')
                                                        #         st.text('Probability is '+ str((probs_nec_list)[i]))
                                                            
                                                        #     elif (pos_ro_nec[i] == "gram+"):
                                                        #         st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-positive ✚ bacteria.')
                                                        #         st.text('Probability is '+ str((probs_poe_list)[i]))
                                                            
                                                        #     elif (pos_ro_nec[i] == "gram-"):
                                                        #         st.success(' 󠀠 󠀠✔️ 󠀠 Resist gram-negative ▬ bacteria.' )
                                                        #         st.text('Probability is '+ str((probs_nec_list)[i]))
                                                            
                                                        #     else:
                                                        #         st.success(' 󠀠 󠀠✔️ 󠀠 Resist other gram of bacteria.')

                                                        # elif anti_or_non[i] == "non antimicrobial":                                                
                                                        #     st.subheader('❌ Your peptide is non antimicrobial peptide.')
                                                        #     st.text('Probability is '+ str((probs_anti_or_non_list)[i]))
                                                        # anti_or_non, pos_ro_nec, probs_anti_or_non_list, probs_nec_list, probs_poe_list = use_model(df_ant_non_normed.iloc[[i]], df_pos_nec_normed.iloc[[i]])
                                                    # feature of peptide --------------------------------------------------------------------------------
                                                    with cooll3:
                                                        
                                                        st.info('🔻 Sequence len: '+ str(len_list[i]))
                                                        st.info('🔻 Hydrophobic: '+ (hydrophobic_list[i]))
                                                        st.info('🔻 Hydrophilic: '+ (hydrophilic_list[i]))
                                                        st.info('🔻 Uncharged: '+ uncharged_list[i])
                                                        #show Positive charge-------------------------------
                                                        st.info('🔻 Positive charge: '+ positive_charge_list[i])

                                                    with cooll5:
                                                            #show Negative charge-------------------------------
                                                        st.info('🔻 Negative charge: '+ Negative_charge_list[i])

                                                            #show Molecular Weight-------------------------------
                                                        st.info('🔻 Molecular Weight: '+ Molecular_Weight_list[i])

                                                            #show Isoelectric Point-------------------------------
                                                        st.info('🔻 Isoelectric Point: '+ pI_list[i])

                                                            #show score hydrophilic-------------------------------
                                                        st.info('🔻 Score hydrophilic: '+ score_hydrophilic_list[i])

                                                            #show Score hydrophobic-------------------------------
                                                        st.info('🔻 Score hydrophobic: '+ Score_hydrophobic_list[i])
                                                    # ("Similarity")
                                                    #     st.info('The similarity of your peptide compare %similarity with another peptide.')
                                                # list_sim_align,  list_iden_align, list_gaps_align = align_sequences(Sequence) list_c_sim, list_matches, list_gaps_al, list_len_al
                                                    with open('style2.css') as f:
                                                        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
                                                        ff01,ff02,ff03,ff04,ff05,ff06,ff07 = st.columns((1,4,0.5,4,0.5,4,1))
                                                        with ff02:
                                                            Ideal_Feature = '<div align="center"><p style="font-family:; color:#1F3D7C; font-size: 18px; background-color: #C2DFFF;">Similarity Comparison</p>'
                                                            st.markdown(Ideal_Feature, unsafe_allow_html=True)
                                                    
                                                            st.write('Defensin peptide is: ' + str(list_c_sim[0]) + '/' + str(list_len_al[0]) + '=' + str(list_sim_align[0])+ '%')    
                                                            st.write('Drosocin peptide is: ' + str(list_c_sim[1]) + '/' + str(list_len_al[1]) + '=' + str(list_sim_align[1])+ '%')
                                                            st.write('Spaetzle peptide is: ' + str(list_c_sim[2]) + '/' + str(list_len_al[2]) + '=' + str(list_sim_align[2])+ '%')
                                                            st.write('B-RAF peptide is: ' + str(list_c_sim[3]) + '/' + str(list_len_al[3]) + '=' + str(list_sim_align[3])+ '%')
                                                            st.write('Hemoglobin peptide is: '+ str(list_c_sim[4]) + '/' + str(list_len_al[4]) + '=' + str(list_sim_align[4])+ '%')
                                                            st.write('Keratin peptide is: '+ str(list_c_sim[5]) + '/' + str(list_len_al[5]) + '=' + str(list_sim_align[5])+ '%')
                                                    # with st.expander("Identity"):
                                                    #     st.info('The identity of your peptide compare identity with another peptide.')
                                                    # list_sim_align,  list_iden_align, list_gaps_align = align_sequences(Sequence) 
                                                        with ff04:
                                                            Ideal_Feature = '<div align="center"><p style="font-family:; color:#1F3D7C; font-size: 18px; background-color: #C2DFFF;">Identity Comparison</p>'
                                                            st.markdown(Ideal_Feature, unsafe_allow_html=True)

                                                            st.write('Defensin peptide is: ' + str(list_matches[0]) + '/' + str(list_len_al[0]) + '=' +  str(list_iden_align[0])+ '%')    
                                                            st.write('Drosocin peptide is: ' + str(list_matches[1]) + '/' + str(list_len_al[1]) + '=' +  str(list_iden_align[1])+ '%')
                                                            st.write('Spaetzle peptide is: ' + str(list_matches[2]) + '/' + str(list_len_al[2]) + '=' +  str(list_iden_align[2])+ '%')
                                                            st.write('B-RAF peptide is: ' + str(list_matches[3]) + '/' + str(list_len_al[3]) + '=' +  str(list_iden_align[3])+ '%')
                                                            st.write('Hemoglobin peptide is: ' + str(list_matches[4]) + '/' + str(list_len_al[4]) + '=' +  str(list_iden_align[4])+ '%')
                                                            st.write('Keratin peptide is: ' + str(list_matches[5]) + '/' + str(list_len_al[5]) + '=' +  str(list_iden_align[5])+ '%')
                                                            
                                                    # with st.expander("Gaps"):
                                                    #     st.info('The gaps of your peptide compare gaps with another peptide.')
                                                    # # list_sim_align,  list_iden_align, list_gaps_align = align_sequences(Sequence)
                                                        with ff06:
                                                            Ideal_Feature = '<div align="center"><p style="font-family:; color:#1F3D7C; font-size: 18px; background-color: #C2DFFF;">Gaps Comparison</p>'
                                                            st.markdown(Ideal_Feature, unsafe_allow_html=True) 

                                                            st.write('Defensin peptide is: ' + str(list_gaps_al[0]) + '/' + str(list_len_al[0]) + '=' +  str(list_gaps_align[0])+ '%')    
                                                            st.write('Drosocin peptide is: ' + str(list_gaps_al[1]) + '/' + str(list_len_al[1]) + '=' +  str(list_gaps_align[1])+ '%')
                                                            st.write('Spaetzle peptide is: ' + str(list_gaps_al[2]) + '/' + str(list_len_al[2]) + '=' +  str(list_gaps_align[2])+ '%')
                                                            st.write('B-RAF peptide is: ' + str(list_gaps_al[3]) + '/' + str(list_len_al[3]) + '=' +  str(list_gaps_align[3])+ '%')
                                                            st.write('Hemoglobin peptide is: ' + str(list_gaps_al[4]) + '/' + str(list_len_al[4]) + '=' +  str(list_gaps_align[4])+ '%')
                                                            st.write('Keratin peptide is: ' + str(list_gaps_al[5]) + '/' + str(list_len_al[5]) + '=' +  str(list_gaps_align[5])+ '%')

                                                    # Graph show amino acid ---------------------------------------------------------------
                                                    c1,c2,c3 = st.columns((4,10,4))
                                                    with c2:
                                                        num_R = df_user_name_seq['Sequence'][i].count('R')
                                                        num_N = df_user_name_seq['Sequence'][i].count('N')
                                                        num_D = df_user_name_seq['Sequence'][i].count('D')
                                                        num_Q = df_user_name_seq['Sequence'][i].count('Q')
                                                        num_E = df_user_name_seq['Sequence'][i].count('E')
                                                        num_K = df_user_name_seq['Sequence'][i].count('K')

                                                        num_H = df_user_name_seq['Sequence'][i].count('H')
                                                        num_S = df_user_name_seq['Sequence'][i].count('S')
                                                        num_T = df_user_name_seq['Sequence'][i].count('T')
                                                        num_A = df_user_name_seq['Sequence'][i].count('A')
                                                        num_C = df_user_name_seq['Sequence'][i].count('C')
                                                        num_F = df_user_name_seq['Sequence'][i].count('F')
                                                        num_G = df_user_name_seq['Sequence'][i].count('G')
                                                        num_I = df_user_name_seq['Sequence'][i].count('I')
                                                        num_L = df_user_name_seq['Sequence'][i].count('L')
                                                        num_M = df_user_name_seq['Sequence'][i].count('M')
                                                        num_P = df_user_name_seq['Sequence'][i].count('P')
                                                        num_V = df_user_name_seq['Sequence'][i].count('V')
                                                        num_W = df_user_name_seq['Sequence'][i].count('W')
                                                        num_Y = df_user_name_seq['Sequence'][i].count('Y')

                                                        df = pd.DataFrame(
                                                        dict(
                                                            hi_ = [num_R, num_N, num_D, num_Q, num_E, num_K, num_H, num_S, num_T, num_A, num_C, num_F, num_G, num_I, num_L, num_M, num_P, num_V, num_W, num_Y],
                                                            bar_labels = ['R','N','D','Q','E', 'K', 'H', 'S', 'T', 'A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y']
                                                        )
                                                        )
                                                        

                                                        df_sorted = df.sort_values('hi_', ascending=False)
                                                                                                
                                                        # def check_freq(x):
                                                        #     freq = {}
                                                        #     for c in set(x):
                                                        #         freq[c] = x.count(c)
                                                        #     return freq

                                                        # freq = check_freq(df_user_name_seq['Sequence'][i])
                                                        # st.write(type(freq))
                                                        # df1 = pd.DataFrame(list(freq.items())).T
                                                        # df1.columns = df1.iloc[0]
                                                                                                    
                                                        # st.dataframe(df1)
                                                        font = {'size': 4}
                                                        # using rc function
                                                        plt.rc('font', **font)
                                                        f, ax = plt.subplots(figsize=(3,1))                                          
                                                        plt.xlabel('amino acid numbers', fontsize=4)
                                                        plt.ylabel('amino acid', fontsize=4)
                                                        plt.title("Amino acid counting diagram", fontsize=5)
                                                        plt.bar('bar_labels', 'hi_', data= df_sorted, color='#1F3D7C')
                                                        
                                                    
                                                        st.pyplot(plt)
                                                    
                                            st.write("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
                                                                                                
                                        # show Dataframe of predict your peptide----------------------------------------------------------------------------------  
                                        df_user_name_seq['Predict_Peptide'] = anti_or_non
                                        df_user_name_seq['Probability_Peptide'] = probs_anti_or_non_list
                                        df_user_name_seq['Resist_Gram'] = pos_ro_nec
                                        df_user_name_seq['Probability_Negative'] = probs_nec_list
                                        df_user_name_seq['Probability_Positive'] = probs_poe_list
                                       

                                        final_data = pd.concat([df_user_name_seq, df_use_in_model], axis=1)
                                        final_data = np.round(final_data, decimals = 5)
                                        if len(final_data) != 0:
                                            html_temp = """
                                                                <div style="background-color:{};height:{};width:{};">
                                                                </div>
                                                                <div style="background-color:#1F3D7C;color:white;padding:10px;border-radius:5px">
                                                                <div id="head" style="background-color:{};padding:1px;border-radius:'15px';">
                                                                </div>
                                                                """
                                            Ideal_dataf = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 40px; background-color: #1F3D7C; border-radius: 5px; text-align:center;">Output</p>'
                                            st.markdown(Ideal_dataf, unsafe_allow_html=True)
                                        
                                            # st.markdown(html_temp.format('white','1px', '50%','#1F3D7C','#FFFFFF'),unsafe_allow_html=True)
                                            st.dataframe(final_data.style.set_properties(**{'background-color': '#C2DFFF',
                                                                                                    'color': 'black',
                                                                                                    'border-color': '#06BBCC',
                                                                                                    }))
                                            with st.container():
                                                left1, left2, centerrr, right1, right2,lasttt = st.columns((8,2,3,5,5,0.3))
                                                with left1:
                                                                    
                                                    def convert_df(final_data):
                                                        return final_data.to_csv(index=False).encode('utf-8')
                                                    csv = convert_df(final_data)
                                                                
                                                    st.download_button(
                                                                "Select to Download .csv file 📥",
                                                                csv,
                                                                "file.csv",
                                                                "text/csv",
                                                                key='download-csv'
                                                                )
                                                    st.write('##')
                                                        
                                                    #color of button --------------------------------------------------------
                                                    m = st.markdown("""
                                                    <style>
                                                    div.stButton > button:first-child {
                                                        background-color: #C2DFFF;
                                                    }
                                                    </style>""", unsafe_allow_html=True)
                                    elif len(df_ant_non_normed) > 50 :
                                        for i in range(len(df_ant_non_normed)):   
                                            anti_or_non, pos_ro_nec, probs_anti_or_non_list, probs_nec_list, probs_poe_list = use_model(df_ant_non_normed.iloc[[i]], df_pos_nec_normed.iloc[[i]])
                                        # show Dataframe of predict your peptide----------------------------------------------------------------------------------  
                                        df_user_name_seq['Predict_Peptide'] = anti_or_non
                                        df_user_name_seq['Probability_Peptide'] = probs_anti_or_non_list
                                        df_user_name_seq['Resist_Gram'] = pos_ro_nec
                                        df_user_name_seq['Probability_Negative'] = probs_nec_list
                                        df_user_name_seq['Probability_Positive'] = probs_poe_list
                                        

                                        final_data = pd.concat([df_user_name_seq, df_use_in_model], axis=1)
                                        final_data = np.round(final_data, decimals = 2)
                                        if len(final_data) != 0:
                                            html_temp = """
                                                                <div style="background-color:{};height:{};width:{};">
                                                                </div>
                                                                <div style="background-color:#1F3D7C;color:white;padding:10px;border-radius:5px">
                                                                <div id="head" style="background-color:{};padding:1px;border-radius:'15px';">
                                                                </div>
                                                                """
                                            Ideal_dataf = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 40px; background-color: #1F3D7C; border-radius: 5px; text-align:center;">Output</p>'
                                            st.markdown(Ideal_dataf, unsafe_allow_html=True)
                                        
                                            # st.markdown(html_temp.format('white','1px', '50%','#1F3D7C','#FFFFFF'),unsafe_allow_html=True)
                                            st.dataframe(final_data.style.set_properties(**{'background-color': '#C2DFFF',
                                                                                                    'color': 'black',
                                                                                                    'border-color': '#06BBCC',
                                                                                                    }))                                       
    
                                        # #download file -------------------------------------------------------------
                                        with st.container():
                                            left1, left2, centerrr, right1, right2,lasttt = st.columns((8,2,3,5,5,0.3))
                                            with left1:
                                                                
                                                def convert_df(final_data):
                                                    return final_data.to_csv(index=False).encode('utf-8')
                                                csv = convert_df(final_data)
                                                            
                                                st.download_button(
                                                            "Select to Download .csv file 📥",
                                                            csv,
                                                            "file.csv",
                                                            "text/csv",
                                                            key='download-csv'
                                                            )
                                                st.write('##')
                                                    
                                                #color of button --------------------------------------------------------
                                                m = st.markdown("""
                                                <style>
                                                div.stButton > button:first-child {
                                                    background-color: #C2DFFF;
                                                }
                                                </style>""", unsafe_allow_html=True)

                                    
                    except:
                        # except Exception as e:
                        # st.error("❌ Failed to load model_anti_or_non")
                        # st.write(f"Details: {str(e)}")   # แสดง error
                        Ideal_err = '<div align="center"><p style="font-sans-serif:; color: white; font-size: 20px; background-color: #F75D59; border-radius: 5px;">Error peptide format!!</p>'
                        st.markdown(Ideal_err, unsafe_allow_html=True)
                        # st.error('Error format your peptide !!')
                        Ideal_forexa = '<div align="left"><p style="font-sans-serif:; color: black; font-size: 15px; background-color: white; border-radius: 5px;">For example, input your peptide FASTA</p>'
                        st.markdown(Ideal_forexa, unsafe_allow_html=True) 
                        # st.write("For example, input your peptide FASTA format")                                  
                        html_temp = """
                                            <div style="background-color:#D1F0FF;padding:1px">
                                            <h8 style="color:black;text-align:left;font-size:80%;"><u>Sample1</u><br>>Sequence_name1<br>DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW<br> 
                                            <br><u>Sample2</u><br>>Sequence_name1<br>DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW<br>>Sequence_name2<br>FPFLLSLIPSAISALKKL </h1>
                                            </div><br>"""
                        st.markdown(html_temp,unsafe_allow_html=True)

                elif cl003.button("❓ how to paste input"):
                    # st.write("For example, input your peptide FASTA format")
                    cffs1,cffs2,cffs3 = st.columns((0.96,12,0.65))
                    with cffs2:
                        Ideal_forexa = '<div align="left"><p style="font-sans-serif:; color: black; font-size: 15px; background-color: white; border-radius: 5px;">For example, input your peptide FASTA</p>'
                        st.markdown(Ideal_forexa, unsafe_allow_html=True)                                  
                        html_temp = """
                                            <div style="background-color:#D1F0FF;padding:1px">
                                            <h8 style="color:black;text-align:left;font-size:80%;"><u>Sample1</u><br>>Sequence_name1<br>DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW<br> 
                                            <br><u>Sample2</u><br>>Sequence_name1<br>DFASCHTNGGICLPNRCPGHMIQIGICFRPRVKCCRSW<br>>Sequence_name2<br>FPFLLSLIPSAISALKKL </h1>
                                            </div><br>"""
                        st.markdown(html_temp,unsafe_allow_html=True)

            
        
