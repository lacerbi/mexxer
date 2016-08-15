# mexxer


### How to make MEX work on your computer

- List of compatible C compilers: http://www.mathworks.com/support/compilers/R2016a/index.html
- Installing a compiler that MATLAB likes can be hard.
  - If you have Windows, you can try with following the instructions [here](https://www.mathworks.com/matlabcentral/answers/101105-how-do-i-install-microsoft-windows-sdk-7-1). Several things can go wrong; see [here](https://www.mathworks.com/matlabcentral/answers/95039-why-does-the-sdk-7-1-installation-fail-with-an-installation-failed-message-on-my-windows-system) for some troubleshooting.
  - As an alternative for Windows users, if you already have a C compiler that you want MATLAB to recognize, download [this script](https://github.com/lacerbi/mexxer/blob/master/mexopts.bat) provided by Bas van Opheusden and follow his instructions: "After downloading this file, move it to *C:\Users\yourusername\AppData\Roaming\MathWorks\MATLAB\yourmatlabversion\*. This folder is hidden but it does exist. Then you'll need to edit the name of the compiler and its path in the mexopts.bat script. 
"

### References

- Bible for MEX files: http://www.mathworks.com/help/matlab/write-cc-mex-files.html
  - C/C++ Matrix Library API: http://www.mathworks.com/help/matlab/cc-matrix-library-api.html
  - Library of example C (and Fortran) files: http://www.mathworks.com/help/matlab/matlab_external/table-of-mex-file-source-code-files.html

### License

This is free, open source software released under the GPL-3 License. See LICENSE.txt for details.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Use this software at your own risk.
