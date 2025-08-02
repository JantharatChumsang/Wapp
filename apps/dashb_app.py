import streamlit as st
import hydralit_components as hc
import codecs
from hydralit import HydraHeadApp
import streamlit.components.v1 as stc 
import base64
import os

class DashbApp(HydraHeadApp):

    def __init__(self, title = 'Dashboard', delay=0, **kwargs):
        self.title = title

    def run(self):
        st.title("Dashboard for data set")
        st.info("Gain a better understanding of the data powering our model. This dashboard reveals trends and relationships used during model training and testing.")

        pdf_path = "apps/Handbook for dashboard.pdf"

        st.markdown(
    """
    <a href="https://shine-tablecloth-195.notion.site/Handbook-for-dashboard-24307ee6772080469424e9e7634df24d?source=copy_link" 
       target="_blank"
       class="external-link-button">
        📘 Open hand-on Dashboard (Notion)
    </a>
    """,
    unsafe_allow_html=True
)
        # Power BI HTML
        def st_webpage(page_html, width=1190, height=600):
            page_file = codecs.open(page_html, 'r')
            page = page_file.read()
            stc.html(page, width=width, height=height, scrolling=False)

        st_webpage('apps/powerBI.html')

        
        
        
        
