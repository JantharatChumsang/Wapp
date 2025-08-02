import streamlit as st
import hydralit_components as hc
import codecs
from hydralit import HydraHeadApp
import streamlit.components.v1 as stc 
import base64

class DashbApp(HydraHeadApp):

    def __init__(self, title = 'Dashboard', delay=0, **kwargs):
        # self.__dict__.update(kwargs)
        self.title = title
        # self.delay = delay

    def run(self):
        st.title("Dashboard for data set")
        st.info("Gain a better understanding of the data powering our model. This dashboard reveals trends and relationships used during model training and testing.")
                #### import html ####

        # Handbook PDF    
        with open("apps/Handbook for dashboard.pdf", "rb") as f:
            pdf_bytes = f.read()
            b64_pdf = base64.b64encode(pdf_bytes).decode()
        
        href = f'''
            <a href="data:application/pdf;base64,{b64_pdf}" target="_blank" class="open-button">
                📄 Open file
            </a>
        '''
        
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
        📄 เปิด Handbook
    </a>
    """,
    unsafe_allow_html=True
)

        
        st.markdown(href, unsafe_allow_html=True)

        
        def st_webpage(page_html,width=1190,height=600):
            page_file = codecs.open(page_html,'r')
            page =page_file.read()
            stc.html(page,width=width, height=height , scrolling = False)
        st_webpage('apps/powerBI.html')
        
        
        
        
