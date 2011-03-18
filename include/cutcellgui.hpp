#ifndef CUTCELLGUI_HPP
#define CUTCELLGUI_HPP

#include <ui_cutcell.h>

namespace cutcell {

namespace gui {

class cutcellgui : public QWidget, private Ui::cutcell
{
    Q_OBJECT

public:
    cutcellgui(QWidget *parent = 0);

public slots:
    void generateSlot();
    void getInFile();
    void getOutFile();
};

} // namespace gui
} // namespace cutcell

#endif // CUTCELLGUI_HPP
