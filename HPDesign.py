import queue
import sys

from PyQt6 import QtGui, QtCore
from PyQt6.QtCore import Qt, QSize, QPoint, QRect
from PyQt6.QtWidgets import QMainWindow, QApplication, QPushButton, QDockWidget, QLabel, QTextEdit, QWidget, QVBoxLayout, QHBoxLayout, QSpinBox, QGroupBox, QDoubleSpinBox, QRadioButton, QComboBox

from BlurWindow.blurWindow import blur

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT

import numpy as np
import scipy.signal as sig
from sympy import *

from control import matlab


import numpy as np

def print_constants(constants : list, var : str):
    NUM = "\n"
    for z in range(1, constants.__len__()+1):
        NUM = NUM + var + str(z) + ' = ' + str(constants[z-1]) + '\n'
        
    return NUM

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)

        self.axes = fig.subplots(2, 1)
        fig.subplots_adjust(left=0.2,
                    bottom=0.2, 
                    right=0.9, 
                    top=0.9, 
                    wspace=1, 
                    hspace=1)

        super(MplCanvas, self).__init__(fig)

class Btn(QPushButton):
    def __init__(self, name: str, x_pos: int, **kwargs):
        super().__init__(**kwargs)

        self.H = queue.Queue(maxsize=1)
        self.name = name
        self.alive = False
        self.setGeometry(QRect(x_pos, 10, 35, 35))

        self.setStyleSheet("""
                            QPushButton
                            {
                                color: white;
                                border: 3px solid rgb(128, 0, 0);
                                image: url(img/""" + name + """.png);
                                border-radius: 10px;
                                background: rgb(0, 0, 0);
                                min-width: 0em;
                            }
                            QPushButton:hover
                            {
                                border: 2px solid rgb(128, 0, 0);
                                border-radius: 10px;
                            }
                            QPushButton:flat
                            { 
                                background: rgb(128, 0, 0);
                            }
                        """)
           

    def enterEvent(self, event: QtGui.QEnterEvent) -> None:
        self.setToolTip(self.name)

    def hitButton(self, pos: QtCore.QPoint) -> bool:
        if self.isDown() and self.name == "Record":
            if not self.alive:
                self.setFlat(True)
                self.alive = True
            else:
                self.alive = False
                self.setFlat(False)

        return True

# Create QDockWidget as an display set
class Doc(QDockWidget):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        self.setFeatures(self.features().DockWidgetMovable)

        self.setStyleSheet("""
                            QDockWidget::title
                            {
                                image: url(img/Praxiscan.jpg);
                                padding-right: 1400px;
                                border-bottom: 2px solid rgba(128, 0, 0, 0.3);
                                padding-bottom: 7px;
                                padding-top: 7px;
                            }
                            QDockWidget::close-button
                            {
                                background: rgba(255, 255, 255, 0);
                                margin: 20px;
                            }
                            
                           """)

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.parentWidget().close()


class MyApp(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self, None)

        self.__press_pos = QPoint()
        self.setGeometry(100, 100, 920, 510)
        # self.setFixedSize(1020, 620)
        
        self.setWindowTitle("EERI 414 - Highpass Design")
        # self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        self.setWindowOpacity(1)
        blur(self.winId())

        doc1 = Doc(parent=self)
        # icon = QtGui.QIcon('./img/ViAxi.ico')
        self.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea , doc1)
        self.setIconSize(QSize(44, 44))
        self.statusBar().showMessage("Chebyshev, Low-Pass IIR, max 0.5 dB loss - Designer")
        self.setStyleSheet("""
                            QStatusBar
                            {
                                background-color : rgba(128, 0, 0, 0.3);
                                border-top: 1px solid rgba(128, 0, 0, 1);
                            }
                            QSizeGrip
                            {
                                image: url(img/Marker.png);
                                width: 16px;
                                height: 16px;
                            }
                           """)
        
        self.Record = Btn("Record", 600 + 177, parent=doc1)
        self.Realization = Btn("Realize Filter", 600 + 230, parent=doc1)
        
        Hlayout = QHBoxLayout()
        Adjustables = QHBoxLayout()
        display = QWidget()
        parWid = QWidget()
        plotWid = QWidget()
        plotType = QGroupBox(plotWid)
        plotType.setTitle("Plot Types:")
        plotLayout = QVBoxLayout()
        parBox = QGroupBox(parWid)
        parBox.setTitle("Design Parameters:")
        parBox.setGeometry(0,0,300,250)
        parameters = QHBoxLayout()
        Left = QVBoxLayout()
        Left.addWidget(QLabel("Fp:"))
        self.Fp = QSpinBox()
        self.Fp.setRange(0,1000000)
        self.Fp.setValue(50000)
        self.Fp.setSingleStep(1000)
        Left.addWidget(self.Fp)
        Left.addWidget(QLabel("Fs:"))
        self.Fs = QSpinBox()
        self.Fs.setRange(0,1000000)
        self.Fs.setValue(80000)
        self.Fs.setSingleStep(1000)
        Left.addWidget(self.Fs)
        Left.addWidget(QLabel("Amax:"))
        self.Amax = QDoubleSpinBox()
        self.Amax.setRange(0.5,0.5)
        self.Amax.setEnabled(False)
        self.Amax.setValue(0.5)
        Left.addWidget(self.Amax)
        Left.addWidget(QLabel("Number of inputs:"))
        self.n = QSpinBox()
        self.n.setRange(0,24)
        self.n.setValue(5)
        Left.addWidget(self.n)
        parameters.addLayout(Left)

        Right = QVBoxLayout()
        Right.addWidget(QLabel("Ft:"))
        self.Ft = QSpinBox()
        self.Ft.setRange(1,1000000)
        self.Ft.setValue(400000)
        self.Ft.setSingleStep(10000)
        Right.addWidget(self.Ft)
        Right.addWidget(QLabel("Order:"))
        self.Order = QSpinBox()
        self.Order.setRange(1,5)
        self.Order.setValue(4)
        Right.addWidget(self.Order)
        Right.addWidget(QLabel("-Amin:"))
        self.Amin = QSpinBox()
        self.Amin.setRange(0,200)
        self.Amin.setValue(40)
        Right.addWidget(self.Amin)
        Right.addWidget(QLabel("Order of Bandwidt:"))
        self.Pow = QSpinBox()
        self.Pow.setRange(0,12)
        self.Pow.setValue(5)
        Right.addWidget(self.Pow)
        parameters.addLayout(Right)

        parBox.setLayout(parameters)

        Adjustables.addWidget(parBox, 80)

        self.Chevy = QRadioButton("Chevy(s)")
        plotLayout.addWidget(self.Chevy)
        self.ALP = QRadioButton("ALP(s)")
        plotLayout.addWidget(self.ALP)
        self.HD = QRadioButton("HD(s)")
        plotLayout.addWidget(self.HD)
        self.GD = QRadioButton("GD(z)")
        plotLayout.addWidget(self.GD)
        self.TD = QRadioButton("TD(YX)")
        plotLayout.addWidget(self.TD)
        self.FDM = QRadioButton("FD|Y|;|X|")
        plotLayout.addWidget(self.FDM)
        self.FDP = QRadioButton("FD/_Y;/_X")
        plotLayout.addWidget(self.FDP)

        plotType.setLayout(plotLayout)

        Adjustables.addWidget(plotType, 20)

        graphLayout = QVBoxLayout()
        self.graph = MplCanvas(self, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar2QT(self.graph, self, coordinates=False)
        graphLayout.addWidget(self.toolbar)
        graphLayout.addWidget(self.graph)

        Hlayout.addLayout(Adjustables, 33)
        Hlayout.addLayout(graphLayout, 66)

        displayLayout = QVBoxLayout()

        displayLayout.addLayout(Hlayout, 75)

        self.info = QTextEdit()
        self.info.setReadOnly(True)
        displayLayout.addWidget(self.info, 25)

        display.setLayout(displayLayout)

        self.setCentralWidget(display)  

        self.Fp.valueChanged.connect(self.updateRequest)
        self.Fs.valueChanged.connect(self.updateRequest)
        self.Amax.valueChanged.connect(self.updateRequest)
        self.Amin.valueChanged.connect(self.updateRequest)
        self.Ft.valueChanged.connect(self.updateRequest)
        self.Order.valueChanged.connect(self.updateRequest)

        self.Chevy.toggled.connect(self.updateRequest)
        self.ALP.toggled.connect(self.updateRequest)
        self.HD.toggled.connect(self.updateRequest)
        self.GD.toggled.connect(self.updateRequest)
        self.TD.toggled.connect(self.updateRequest)
        self.FDM.toggled.connect(self.updateRequest)
        self.FDP.toggled.connect(self.updateRequest)

        self.Realization.pressed.connect(self.realize)

        self.Chevy.setChecked(True)
        self.Record.alive = True
        self.Record.setFlat(True)
        self.updateRequest()

    def realize(self):
        w_p = (2*np.pi*self.Fp.value())/self.Ft.value()
        w_s = (2*np.pi*self.Fs.value())/self.Ft.value()

        Omega_p = np.tan(w_p/2)
        Omega_s = np.tan(w_s/2)
        BIG_omega = Omega_s/Omega_p
        
        if self.Order.value() == 1:
            myZeros = [4.0811]
            myPoles = [1, 4.10811]
        if self.Order.value() == 2:
            myZeros = [2.05405]
            myPoles = [1, 1.79668, 2.11403]
        if self.Order.value() == 3:
            myZeros = [1.022702]
            myPoles = np.polymul([1, 0.76722, 1.33863], [1, 0.76722])
        if self.Order.value() == 4:
            myZeros = [0.51352]
            myPoles = np.polymul([1, 0.42504, 1.16195], [1, 1.02613, 0.45485])
        elif self.Order.value() == 5:
            myZeros = [0.25676]
            myPoles = np.polymul([1, 0.27005, 1.09543], [1, 0.70700, 0.53642])
            myPoles = np.polymul(myPoles,[1,0.43695])

        p = []

        for i in np.arange(self.Order.value(), -1, -1):
            p.append(myPoles[i]*np.power(Omega_p,self.Order.value()-i)/myPoles[self.Order.value()])
        
        p.append(myZeros[0]/myPoles[self.Order.value()])

        N = [2*self.Ft.value(),-2*self.Ft.value()]			
        D = [1, 1]			
        HPDEN = [0]

        for i in range(p.__len__() - 1):
            DEN = [1]
            for _ in np.arange(i, self.Order.value(), 1):  
                DEN = np.polymul(D,DEN)

            for _ in np.arange(i):
                DEN = np.polymul(N,DEN)

            if i == p.__len__() - 2:
                D4 = np.polymul(p[p.__len__() - 1],DEN)
            
            DEN = np.polymul(p[i],DEN) # use p1 = p[self.Order.value()] for D^(self.Order.value())
            HPDEN = np.polyadd(HPDEN,DEN) 

        D4 = D4/HPDEN[0]
        HPDEN = np.flip(HPDEN/HPDEN[0])  

        A = D4      # [p4, p3, p2, p1] all 'a' values (not all pass)
        N = self.Order.value()
        for n in range(1,N + 1):
            A[n] = A[n] - A[0]*HPDEN[N-n+1] # High to low [a4, a3, a2...]

        cnt = 1 # A will always have A[0] and A[1] since 1 <= N <= 5
        DSET = HPDEN # low to high [1, d1, d2, d3, d4]
        S = []
        S = D4
        d = HPDEN
        k = []  # High to low [k4, k3, k2...]
        k.append(D4[0])  # Alpass filter d0
        for n in np.arange(N,1,-1): 
            d = [] 
            for i in np.arange(1,n): # 1 2 3
                d.append((DSET[i]-DSET[0]*DSET[n-i])/(1-DSET[0]**2)) # High to low [d3', d2', d1']
            DSET = d            #          >--v--< k is always first highest coefficient
            k.append(DSET[0])  # High to low [k, d3, d2...]
            S = np.polyadd(S,np.polymul(DSET, A[N-n]))
            cnt += 1
            for r in range(cnt,N+1):
                A[r] = A[r] - A[cnt-1]*DSET[r-cnt]
        
        CaseInfo = "\n_______________________________________________________________________________________________________________________________________________\n"
        CaseInfo = CaseInfo + "\nGrey-Markel realization\n" + print_constants(A,'a') + print_constants(np.flip(k),'k')

        self.info.setText('|w_p: ' + str(w_p) + ' rad/s\t|w_s: ' + str(w_s) +' rad/s\t|\n|Omega_p: '
            + str(Omega_p) + ' Units\t|Omega_s: ' + str(Omega_s) + ' Units\t|BIG_omega: '+ str(BIG_omega) +' Unitless'
            + CaseInfo)
            

    def updateRequest(self):
        if self.Record.alive == True:
            CaseInfo = "\n_______________________________________________________________________________________________________________________________________________\n"

            self.graph.axes[0].clear()
            self.graph.axes[1].clear()

            Omega_p = 2*np.pi*self.Fp.value()
            Omega_s = 2*np.pi*self.Fs.value()

            w_p = Omega_p/self.Ft.value()
            w_s = Omega_s/self.Ft.value()

            BIG_omega = Omega_s/Omega_p

            epsalon = (10**(0.1*0.25)-1)**0.5

            if self.Order.value() == 1:
                DEN = np.poly1d([0, 2.86278])
                NUM = np.poly1d([1, 2.86278])
            if self.Order.value() == 2:
                DEN = np.poly1d([0, 0, 1.43138])
                NUM = np.poly1d([1, 1.42562, 1.51620])
            if self.Order.value() == 3:
                DEN = np.poly1d([0, 0, 0, 0.71570])
                NUM = np.poly1d(np.polymul([1, 0.62646, 1.14245], [1, 0.62646]))
            if self.Order.value() == 4:
                DEN = np.poly1d([0, 0, 0, 0, 0.35785])
                NUM = np.poly1d(np.polymul([1, 0.35071, 1.06352], [1, 0.84668, 0.356412]))
            elif self.Order.value() == 5:
                DEN = np.poly1d([0, 0, 0, 0, 0, 0.25676])
                NUM = np.poly1d(np.polymul([1, 0.22393, 1.03578], [1, 0.58625, 0.47677]))
                NUM = np.poly1d(np.polymul(NUM,[1,0.362332]))
            
            
            if self.Chevy.isChecked():
            # xxxxxxxxxxxxxxxxx Analog Low-pass xxxxxxxxxxxxxxxxx
                C = sig.lti(NUM, DEN)
                # print('lowpass:', Lowpass)
                w, m, p = sig.bode(C)

                if self.Amin.value() > np.min(m[w>BIG_omega]):
                    self.Amin.setStyleSheet( "QSpinBox { background-color: red; }")
                else:
                    self.Amin.setStyleSheet( "QSpinBox { background-color: white; }")

                self.graph.axes[0].plot(w, m, color = "indigo")
                self.graph.axes[0].set_xscale("log")
                self.graph.axes[0].set_title("Analog Low-pass Magnitude")
                self.graph.axes[0].set_ylabel("dB")
                self.graph.axes[0].set_xlabel("Frequency(rad/s)")
                self.graph.axes[0].grid(1)

                self.graph.axes[1].plot(w, p, color = "red")
                self.graph.axes[1].set_xscale("log")
                self.graph.axes[1].set_title("Analog Low-pass phase")
                self.graph.axes[1].set_xlabel("Frequency(rad/s)")
                self.graph.axes[1].set_ylabel("Degrees")
                self.graph.axes[1].grid(1)

                CaseInfo = CaseInfo + '\nAmin = ' + str(np.min(m[w>BIG_omega])) + ' dB'
                G = matlab.tf(NUM.coeffs, DEN.coeffs)
                CaseInfo = CaseInfo + "\nStep 1.) Normalized Analog Filter=\n" + str(G)

            # xxxxxxxxxxxxxxxxx Analog Low-pass phase xxxxxxxxxxxxxxxxx
            elif self.ALP.isChecked():
            # xxxxxxxxxxxxxxxxx Analog High-pass phase xxxxxxxxxxxxxxxxx

                designNum = np.poly1d([1/Omega_p, 0]) # basically, s = (^s)/Omeg_p to achieve cheby's low pass
                designDEN = np.poly1d([1])

                # all chevyshev numinators are of size 1 in a list
                
                HaNUM, HaDEN = substitude(designNum,designDEN,NUM,z)

                ALP_NUM = HaDEN
                ALP_DEN = np.polyadd(HaDEN,np.polymul(HaNUM,epsalon))
            
                ALP = sig.lti(ALP_NUM, ALP_DEN)
                w, m, p = sig.bode(ALP)
                self.graph.axes[0].plot(w/(2*np.pi), m, color = "indigo")
                self.graph.axes[0].set_xscale("log")
                self.graph.axes[0].set_title("Analog High-pass Magnitude")
                self.graph.axes[0].set_ylabel("dB")
                self.graph.axes[0].set_xlabel("Frequency(Hz)")
                self.graph.axes[0].grid(1)

                self.graph.axes[1].plot(w/(2*np.pi), p, color = "red")
                self.graph.axes[1].set_xscale("log")
                self.graph.axes[1].set_title("Analog High-pass phase")
                self.graph.axes[1].set_xlabel("Frequency(Hz)")
                self.graph.axes[1].set_ylabel("Degrees")
                self.graph.axes[1].grid(1)

                G = matlab.tf(list(ALP_NUM), list(ALP_DEN))
                CaseInfo = CaseInfo + "\nStep 2.) Denormalized Analog Filter=\n" + str(G)
            # xxxxxxxxxxxxxxxxx Analog High-pass phase xxxxxxxxxxxxxxxxx
            elif self.HD.isChecked():
            # xxxxxxxxxxxxxxxxx Digital High-pass phase xxxxxxxxxxxxxxxxx

                designNum = np.poly1d([1/Omega_p, 0]) # basically, s = (^s)/Omeg_p to achieve cheby's low pass
                designDEN = np.poly1d([1])

                # all chevyshev numinators are of size 1 in a list
                
                HaNUM, HaDEN = substitude(designNum,designDEN,NUM,DEN)

                ALP_DEN = np.polyadd(HaDEN,np.polymul(HaNUM,epsalon))
                ALP_NUM = HaDEN   

                fTransNUM = np.poly1d([Omega_s*Omega_p])    # basically, s = Omeg_s/(^s) to transform to high pass
                fTransDEN = np.poly1d([1, 0])

                HDNUM, HDDEN = substitude(fTransNUM,fTransDEN,ALP_NUM,ALP_DEN)

                Hd = sig.TransferFunction(HDNUM, HDDEN)
                
                w, m, p = sig.bode(Hd)
                self.graph.axes[0].plot(w/(2*np.pi), m, color = "indigo")
                self.graph.axes[0].set_xscale("log")
                self.graph.axes[0].set_title("Digital High-pass Magnitude")
                self.graph.axes[0].set_ylabel("dB")
                self.graph.axes[0].set_xlabel("Frequency(Hz)")
                self.graph.axes[0].grid(1)

                self.graph.axes[1].plot(w/(2*np.pi), p, color = "red")
                self.graph.axes[1].set_xscale("log")
                self.graph.axes[1].set_title("Digital High-pass phase")
                self.graph.axes[1].set_xlabel("Frequency(Hz)")
                self.graph.axes[1].set_ylabel("Degrees")
                self.graph.axes[1].grid(1)

                G = matlab.tf(list(HDNUM), list(HDDEN))
                CaseInfo = CaseInfo + "\nStep 3.) Digital Filter=\n" + str(G)

            # xxxxxxxxxxxxxxxxx Digital High-pass phase xxxxxxxxxxxxxxxxx
            elif self.GD.isChecked():

                designNum = np.poly1d([1/Omega_p, 0]) # basically, s = (^s)/Omeg_p to achieve cheby's low pass
                designDEN = np.poly1d([1])

                # all chevyshev numinators are of size 1 in a list
                
                HaNUM, HaDEN = substitude(designNum,designDEN,NUM,DEN)

                ALP_DEN = np.polyadd(HaDEN,np.polymul(HaNUM,epsalon))
                ALP_NUM = HaDEN   

                fTransNUM = np.poly1d([Omega_s*Omega_p]) 
                fTransDEN = np.poly1d([1, 0])

                HDNUM, HDDEN = substitude(fTransNUM,fTransDEN,ALP_NUM,ALP_DEN)   

                z,p = sig.bilinear(HDNUM, HDDEN,fs=self.Ft.value())
                GD = sig.TransferFunction(z,p,dt=1/self.Ft.value())
                
                w, m, p = sig.dbode(GD)
                self.graph.axes[0].plot(w/(2*np.pi), m, color = "indigo")
                self.graph.axes[0].set_xscale("log")
                self.graph.axes[0].set_title("Digital High-pass Magnitude")
                self.graph.axes[0].set_ylabel("dB")
                self.graph.axes[0].set_xlabel("Frequency(Hz)")
                self.graph.axes[0].grid(1)

                self.graph.axes[1].plot(w/(2*np.pi), p, color = "red")
                self.graph.axes[1].set_xscale("log")
                self.graph.axes[1].set_title("Digital High-pass phase")
                self.graph.axes[1].set_xlabel("Frequency(Hz)")
                self.graph.axes[1].set_ylabel("Degrees")
                self.graph.axes[1].grid(1)

                G = matlab.tf(list(z), list(p),dt=1/self.Ft.value())
                CaseInfo = CaseInfo + "\nStep 3.) Digital Filter=\n" + str(G)

            elif self.TD.isChecked():

                pass

            elif self.FDM.isChecked():

                pass

            elif self.FDP.isChecked():

                pass

                

            self.info.setText('|w_p: ' + str(w_p) + ' rad/s\t|w_s: ' + str(w_s) +' rad/s\t|\n|Omega_p: '
            + str(Omega_p) + ' Units\t|Omega_s: ' + str(Omega_s) + ' Units\t|BIG_omega: '+ str(BIG_omega) +' Unitless'
            + CaseInfo)

            self.graph.draw()
            
                
    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.__press_pos = event.pos()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.__press_pos = QPoint()

    def mouseMoveEvent(self, event):
        if not self.__press_pos.isNull():
            self.move(self.pos() + (event.pos() - self.__press_pos))

def equalOrder(NUM : list, DEN:list, polyReturnIfEqual : list) -> list:
    NUM = list(NUM)
    DEN = list(DEN)
    polyReturnIfEqual = list(polyReturnIfEqual)

    if np.size(NUM) < np.size(DEN):
        for _ in range(np.size(NUM), np.size(DEN)):
            NUM.append(0)
        return list(np.flip(NUM))
    elif np.size(NUM) > np.size(DEN):
        for _ in range(np.size(DEN), np.size(NUM)):
            DEN.append(0)
        return list(np.flip(DEN))
    
    return list(np.flip(polyReturnIfEqual))
    

def substitude(sub_num : list,sub_den : list,num : list, den : list):
    sub_num = list(sub_num)
    sub_den = list(sub_den)
    num = list(num)
    den = list(den)

    NUM = [0]
    DEN = [0]
    if den.__len__() > num.__len__():
        N = den.__len__()
        num = equalOrder(num,den,num)
    else:
        N = num.__len__()
        den = equalOrder(num,den,den)
    
    for k in range(N):
        Nn = [1]
        Dn = [1]
        for _ in range(N-k - 1):
            Nn = np.polymul(sub_num,Nn)
            Dn = np.polymul(sub_num,Dn)

        Nd = [1]
        Dd = [1]
        for _ in range(k):
            Nd = np.polymul(sub_den,Nd)
            Dd = np.polymul(sub_den,Dd)    
        
        TermNum = np.polymul(Nn,Nd)
        TermDen = np.polymul(Dn,Dd)
        if k == N - 1:
            NUM = np.polyadd(NUM, np.polymul(TermNum,num[k]))
            DEN = np.polyadd(DEN, np.polymul(TermDen,den[k]))
            NUM = np.polymul(NUM,TermDen)
            DEN = np.polymul(DEN,TermNum)
        else:
            NUM = np.polyadd(NUM, np.polymul(TermNum,num[k]))
            DEN = np.polyadd(DEN, np.polymul(TermDen,den[k]))

        
    return NUM,DEN
            
if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle("Windows")
    mywindow = MyApp()
    mywindow.show()
    app.exec()