# .............................................................................................................................
from hydralit import HydraApp
import hydralit_components as hc
import apps
from apps import (
    home,
    predict,
    dashb_app,
    intro_app_new,
    member_app,
    ContactUsApp
)
import apps.HowToApp as how_to_app
import streamlit as st
from streamlit_option_menu import option_menu
from streamlit_lottie import st_lottie
import streamlit.components.v1 as stc
import json
import pandas as pd
import numpy as np
# -------------------------------------------------------------------------------------------------------------------------------------

#Only need to set these here as we are add controls outside of Hydralit, to customise a run Hydralit!
st.set_page_config(page_title='antimicrobial peptide',page_icon=":pill:",layout='wide',initial_sidebar_state='auto',)

if __name__ == '__main__':
    with open('style2.css') as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
        #---ONLY HERE TO SHOW OPTIONS WITH HYDRALIT - NOT REQUIRED, use Hydralit constructor parameters.
        st.write("##")
        hide_st = ('Hide Streamlit Markers',True)
        over_theme = {'txc_inactive': '#FFFFFF','menu_background':'#1F3D7C','txc_active':'#F36C23','option_active':'#FFFFFF'}
        #this is the host application, we add children to it and that's it!
        app = HydraApp(
            title='Secure Hydralit Data Explorer',
            favicon=":pill:",
            hide_streamlit_markers=hide_st,
            #add a nice banner, this banner has been defined as 5 sections with spacing defined by the banner_spacing array below.
            use_banner_images=["resources/pig0.png",None,{'header':"<h1 style='text-align:center;padding: 0px 0px;color:F36C23;font-size:175%;'>WAAPP: Web Application for Antimicrobial Peptide Prediction</h1>"},None,"./resources/lock.png"], 
            banner_spacing=[6,15,60,15,4.5],
            navbar_theme=over_theme
        )

        with open('style2.css') as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
        #Home button will be in the middle of the nav list now
         

        app.add_app("Home", icon="🏠", app=home.HomeApp(title='Home'), is_home=True)
        app.add_app("Predict", icon="🔬", app=predict.PredictApp(title='Predict'))
        app.add_app("Predict your peptide", icon="🔍", app=predict.PredictApp(title="Predict"))

        # add all application classes
        app.add_app("How to use web application", icon="❓", app=how_to_app.HowToApp(title="How to use web application"))
        app.add_app("Dashboard", icon="📊", app=dashb_app.DashbApp(title="Dashboard"))
        app.add_app("Intro", icon="🥂", app=intro_app_new.IntroApp(title="About us"))
        app.add_app("Member", icon="🧑‍💼", app=member_app.MemberApp(title="Member"))
        app.add_app("Contact us", icon="📧", app=ContactUsApp.ContactUsApp(title="Contact us"))

        #check access
        username = app.check_access()

        # def st_webpage(page_html,width=1370,height=1550):
        #     page_file = codecs.open(page_html,'r')
        #     page =page_file.read()
        #     stc.html(page,width=width, height=height , scrolling = False)

        # If the menu is cluttered, just rearrange it into sections!
        # completely optional
        if username:
            complex_nav = {
                'Home': ['Home'],
                'Predict your peptide': ["Predict your peptide"],
                'How to WebApp': ["How to use web application"],
                'Dashboard': ['Dashboard'],
                '🕮 About us': ['Intro',"Member"],
                'Contact us': ['Contact us']
            }
    
        else:
            complex_nav = {
                'Home': ['Home'],
            }

        #and finally just the entire app and all the children.
        app.run(complex_nav)

