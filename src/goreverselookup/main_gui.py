# The main GUI class of GOReverseLookup. Using PyQt6.

import PyQt6.QtWidgets as qtw
import PyQt6.QtGui as qtg
import sys

class MainWindow(qtw.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GOReverseLookup")
        self.setLayout(qtw.QVBoxLayout())
        
        hw_label = qtw.QLabel("Hello world, what's your name?")
        hw_label.setFont(qtg.QFont('Helvetica', 18))
        self.layout().addWidget(hw_label)
        
        # entry box
        entry = qtw.QLineEdit()
        entry.setObjectName("name_field")
        entry.setText("")
        self.layout().addWidget(entry)
        
        # button
        button = qtw.QPushButton("Press me", clicked = lambda: press_it())
        self.layout().addWidget(button)

        self.show() # shows the window
        
        def press_it():
            hw_label.setText(f"Hello, {entry.text()}")
            entry.setText("")


app = qtw.QApplication([])
main_window = MainWindow()
# run the app
app.exec()


    