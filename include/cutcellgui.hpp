#ifndef CUTCELLGUI_HPP
#define CUTCELLGUI_HPP

#include <ui_cutcell.h>

class cutcellgui : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

public:
    cutcellgui(QWidget *parent = 0);

public slots:
    void generateSlot();
    void getInFile();
    void getOutFile();
    void cubicGridSlot();
};

#endif // CUTCELLGUI_HPP
