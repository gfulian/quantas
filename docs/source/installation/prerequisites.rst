.. _prerequisites:

=============
Prerequisites
=============

In order to install and run Quantas, the following software is required:

- python_ >=3.5 (programming language used to develop Quantas)
- python3-pip_ (Python3 package manager)

.. _Python: http://www.python.org/
.. _python3-pip: https://packaging.python.org/tutorials/installing-packages/#requirements-for-installing-packages

Depending on your set up, the following few optional dependencies may be useful:

- virtualenv (software to create a virtual python environment to install Quantas in)
- virtualenvwrapper (a wrapper package for virtualenv, used to easily handle virtual
  environments)

Supported operating systems
===========================
Quantas has been developed in order to be used on several possible operating systems:

- Unix (tested on Debian, Ubuntu and Raspbian)
- Mac OS
- Windows (tested on 7, 8.1 and 10)
- Windows Subsystem for Linux (WSL, tested on 
  `Ubuntu <https://www.microsoft.com/en-gb/p/ubuntu/9nblggh4msv6?source=lp&activetab=pivot:overviewtab>`_)
  
It is expected that Quantas can also run on:

- Older and newer Debian/Ubuntu versions
- Other Linux distributions

.. _Debian: https://www.debian.org/index.it.html
.. _Raspbian: https://www.raspberrypi.org/downloads/raspbian/


Debian/Ubuntu
=============

You can use the ``apt`` package manager to install the prerequisites on Debian and derived 
distributions. The following will install the basic python requirements on your system::

  sudo apt-get install python3-dev python3-pip virtualenv


Mac OS (Homebrew)
=================

Installation of Quantas has been tested using the `Homebrew <http://brew.sh/index_de.html>`_
package manager. If Homebrew is not present in your system, you can install it using the 
following command::

  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
  
After the installation of Homebrew, you can install the prerequisites of Quantas with the following command::

  brew install python


Windows Subsystem for Linux (WSL)
=================================

If you are under Windows 10, you can use WSL to install a Linux distribution of your choice and follow the guide for the Debian/Ubuntu installation of the prerequisites for Quantas.
