#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QFileInfo>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_zip_clicked();

    void on_pushButton_input_zip_clicked();

    void on_pushButton_output_zip_clicked();

    void on_pushButton_input_unzip_clicked();

    void on_pushButton_output_unzip_clicked();

    void on_pushButton_unzip_clicked();


private:
    Ui::MainWindow *ui;

private:
    QFileInfo  _input_zip_file;
    QFileInfo  _input_unzip_file;
    QFileInfo  _output_zip_file;
    QFileInfo  _output_unzip_file;
};
#endif // MAINWINDOW_H
