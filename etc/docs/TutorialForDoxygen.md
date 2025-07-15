# Tutorial for using Doxygen with Fortran Coding
#预警中心 #Doxygen #开发规范

Created by ([Zilong QIN](zilong.qin@gmail.com)), 2020/5/9, @GBA-MWF, Shenzhen

## Doxygen? 
 An Automatic Documentation Tool.
Write DXG style comments in code -> Generating documentation automatically  

## Why Doxygen?
* Synchronous mode, always keep the documentation updated
* **Fortran** -> major coding language
* Best choice - a cross platform tool [http://fortranwiki.org/fortran/show/Automatic+documentation](http://fortranwiki.org/fortran/show/Automatic+documentation)

## Where and how to get?
Download [ http://www.doxygen.nl/download.html](http://www.doxygen.nl/download.html) 
Official manual [http://www.doxygen.nl/manual/index.html](http://www.doxygen.nl/manual/index.html) 
Source code [https://github.com/doxygen/doxygen](https://github.com/doxygen/doxygen)
GraphViz (Required for UML flowchart) [http://www.graphviz.org/](http://www.graphviz.org/)
Homebrew (Package manager)[https://brew.sh/](https://brew.sh/)

### Installation
apt-get:
```bash
apt-get install doxygen
apt-get install graphviz
```

Homebrew:
```bash
brew install doxygen
brew install graphviz
```

## Template for beginner
The template can be achieved through git:
```
git clone https://e.coding.net/opensimi/DTFC.git
```

### Basic rules
DXG comment starts with !> and continues with !!, ends with a line without the previous two.
Tags are identified with @.   @brief, @author, etc.
Markdown style is supported  [http://www.doxygen.nl/manual/markdown.html](http://www.doxygen.nl/manual/markdown.html).

* For a **program** 
	* Universal info title (title, version, history, etc. )
	* DXG style code blocks (briefs, copyrights, etc. )
	* @test 
```
!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           : 
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/12/29, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved. 
! @note 
! @warning 
! @attention 
```

* For a **module**
	* 	* Universal info title (title, version, history, etc. )
	* DXG style code blocks (briefs, copyrights, etc. 

```
!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           : 
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/12/29, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved. 
!! @note 
!! @warning 
!! @attention 
```

* For a **type**
	* Doxygen style code blocks
	* Explanation  for each member

```
!> @brief
!! An class for diffIcos_t
!! @see       model_states
!! @note      a little side note
!! @warning   warning information
!! @attention key information
TYPE, EXTENDS(differentiator_t) :: diffIcos_t
CONTAINS
   PROCEDURE :: derivative !< abstract function for derivative of differentiator
END TYPE diffIcos_t
```

* For a **interface**
```
    !> @brief
    !! The inferface for derivative function in differentiator_t class
    !! @see       model_states
    !! @note      a little side note
    !! @warning   warning information
    !! @attention key information
    ABSTRACT INTERFACE
        FUNCTION diff(this, f, dx)
            IMPORT differentiator
            CLASS(differentiator), INTENT(IN) :: this !< parent reference
            REAL, INTENT(IN) :: f(:)                  !< input function array
            REAL, INTENT(IN) :: dx                    !< d value
            REAL, DIMENSION(SIZE(f)) :: diff          !< size of diff
        END FUNCTION diff
    END INTERFACE
```

* For a **subroutine**

```
!> @brief
!! Build the restriction matrix for the aggregation 1
!! @see       intrestbuild2 -> this can be a hyperlink by matching the function name.
!! @note      a little side note
!! @warning   warning information, optional
!! @attention key information, optional
SUBROUTINE intrest_build1(A, aggr, Restrict, A_ghost)
    IMPLICIT NONE
    TYPE(SpMtx), INTENT(IN) :: A         !< fine level matrix
    TYPE(Aggrs), INTENT(IN) :: aggr      !< aggregation matrix
    TYPE(SpMtx), INTENT(OUT) :: Restrict !< restriction matrix
    !...
END SUBROUTINE
```

* For a **function**
	* 	DXG style code blocks
	*  Explanation  for each input/output argument

```
!> @brief
!! Build the restriction matrix for the aggregation 1
!! @see       intrestbuild2 -> this can be a hyperlink by matching the function name.
!! @note      a little side note
!! @warning   warning information, optional
!! @attention key information, optional
SUBROUTINE intrest_build1(A, aggr, Restrict, A_ghost)
    IMPLICIT NONE
    TYPE(SpMtx), INTENT(IN) :: A         !< fine level matrix
    TYPE(Aggrs), INTENT(IN) :: aggr      !< aggregation matrix
    TYPE(SpMtx), INTENT(OUT) :: Restrict !< restriction matrix
    !...
END SUBROUTINE
```

## Generating documentation with DXG
Basic software routine.
**Keynote**: 
EXTRACT_*
SORT_*
GENERATE_*
HIDE_UNDOC_*                 : Hide the functions and classed without doxygen comments.
SOURCE_BROWSER               : Enable the reference to source files
REFERENCED_BY_RELATION
REFERENCES_RELATION
GENERATE_TREEVIEW            : Add a sider on the left for navigation
HTML_DYNAMIC_*               : Not so useful 

GRAPH
	HAVE_DOT
	DOT_PATH

Remember to save the configuration.

## Where is the Main page, Get started?
Put the README.MD as the main page.

INPUT -> USE_MDFILE_AS_MAINPAGE
README.MD

## Examples
OpenCV [https://docs.opencv.org/4.3.0/](https://docs.opencv.org/4.3.0/)
Biogears [Biogears](https://biogearsengine.com/documentation/index.html)










