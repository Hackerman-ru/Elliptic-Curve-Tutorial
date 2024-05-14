# ����������� �� ���������� ������������ �� ������������� ������
## ��������
������ ������ �������� ������������ � **��** ������ � ����� ������������� ����� submodule.
���� ������ ������������ �������� ���, �� ������ ���������� ��� � ��� ������.
## ������
### ������������
- [Visual Studio 2022](https://visualstudio.microsoft.com/downloads/)
### ��������� ���������
1. ���������� Visual Studio 2022. �� ����� ��������� ��� ����� � Visual Studio Installer ��������:
	- ���������� ������������ ���������� �� C++
	- ������� ������ ��� Google Test
2. ��� ��������� ������������ � ������� "Options -> Test Adapter for Google Test -> Parallelization" ���������� 10 ������ � �������� ������������ ���������� ������
### ������ ��� ������������
1. ����������� ����������� ������ � submodules:
```
git clone --recursive https://github.com/Hackerman-ru/Elliptic-Curve-Tutorial.git
```
2. �������� ����� Visual Studio ������� � ����� `Elliptic-Curve-Tutorial\Elliptic-Curve-Tutorial.sln`
3. �������� �� ������� ������ Release ��� Debug ������������ ������
4. ������� �� ������� ������: Build -> Build Solution
5. ������� �� ������� ������: Test -> Test Explorer
6. ��������� �����, ����� �� "Run All Tests in View"
7. ���� ������ �������������� ����������, ��������� ������� ���������� �� boost ������ ����������, � ����� src/uint.h ���������������� ������
`#define ECG_USE_BOOST` � ������������ �������
## ����������� ����������
<img alt="UML diagram" src="ClassDiagram.png" />