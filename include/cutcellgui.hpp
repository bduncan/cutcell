#ifndef CUTCELLGUI_HPP
#define CUTCELLGUI_HPP

#include <ui_cutcell.h>
#include <QtCore/QThread>

class GenerateThread : public QThread
{
    Q_OBJECT
public:
    GenerateThread(QWidget *parent = 0) : QThread(parent) {};
    void run();
};

class cutcellgui : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

    friend class GenerateThread; // Allow this thread to access the widgets.

public:
    cutcellgui(QWidget *parent = 0);
    QLabel *statusBarLabel;

public slots:
    void generateSlot();
    void getInFile();
    void getOutFile();
    void cubicGridSlot();
    void cubicCellSlot();
    void cubicOriginSlot();

private:
    GenerateThread *thread;
};

#endif // CUTCELLGUI_HPP
