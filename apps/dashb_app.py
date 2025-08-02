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
        <style>
        .open-button {
            display: inline-block;
            padding: 10px 24px;
            font-size: 16px;
            color: #333;
            border: 2px solid #ccc;
            border-radius: 999px;
            text-decoration: none;
            background-color: #fff;
            font-weight: 500;
            transition: background-color 0.3s ease, color 0.3s ease;
        }
        .open-button:hover {
            background-color: #28a745;
            color: white;
            border: 2px solid #28a745;
        }
        </style>
    
        <a href="apps/Handbook%20for%20dashboard.pdf" target="_blank" class="open-button">
            📄 เปิดดู PDF คู่มือ Dashboard
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

        
        
        
        
