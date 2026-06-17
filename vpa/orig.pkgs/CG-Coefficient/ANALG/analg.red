%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                    %%
%%                                                                    %%
%%     PROGRAM-PACKAGE FOR ANGULAR-MOMENTUM ALGEBRA: ANALG            %%
%%                                                                    %%
%%                                         Coded by Fumihiro Koike    %%
%%                                                                    %%
%%                                         April 1992, Version 1.1    %%
%%                                                                    %%
%%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% The full path name of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter00.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter00.r,v 1.1 92/04/16 20:37:35 koikef Locked $'/
%
%
% Contents 
%
%   Prologue   Define basic operators and mode setters
%   Chapter 1. Collection of auxiliary operators and procedures
%   Chapter 2. Clebsch-Gordan coefficients and related quantities
%   Chapter 3. Useful quantities related to the C. G. coefficients
%   Chapter 4. Racah coefficients and Wigner's 6j-Symbols
%   Chapter 5. 9j-symbols or X-coefficients
%   Appendix   Diagnostics routines
%   Epilogue   Initialize global variables
%
% Notes
%
%   (1) The user should not refer any symbols with under bars: "_"
%       directly from the command lines. Otherwise, unexpected
%       change of the system variables could take place.
%       Especially, any substitutions into these symbols
%       should be strictly avoided.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% Prologue.  Define basic operators and mode setters.                  %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Operators :
%
define is = = $
define be = = $
define isnot = neq $
define benot = neq $
define takes = := $
%
%
% Mode setters :
%
%
% If mode!_id!_ is
%   activate: perform formula reduction,
% deactivate: avoid   formula reduction.
%
procedure evaluation ( mode!_id!_ ) $
begin $
  if mode!_id!_ is deactivate then symbolic!_mode!_ takes   active $
  if mode!_id!_ is   activate then symbolic!_mode!_ takes inactive $
  if mode!_id!_ is       show then 
                         write "Currently ", symbolic!_mode!_, "." $
  if mode!_id!_ isnot deactivate and mode!_id!_ isnot activate then
     write "Usage: evaluation [ activate, deactivate, show ]" $
  return symbolic!_mode!_ $
  end $
%
%
% If mode!_id!_ is
% Rose  : transform Wigner's notations to Rose's notations,
% Wigner: transform Rose's notations to Wigner's notations,
% Asitis: avoid transformation.
%
% Note: This feature is effective only in the symbolic mode,
%       i.e., only when a function call
%             "evaluation deactivate"
%       is preceded.
%
procedure notation ( mode!_id!_ ) $
begin scalar wrk0 $
  if mode!_id!_ is Rose then do
    begin $
      c!_to!_w3j!_     takes inactive $
      w3j!_to!_c!_     takes   active $
      w!_to!_w6j!_     takes inactive $
      w6j!_to!_w!_     takes   active $
      wrk0 takes mode!_id!_ $
      end $
  if mode!_id!_ is Wigner then do
    begin $
      c!_to!_w3j!_     takes   active $
      w3j!_to!_c!_     takes inactive $
      w!_to!_w6j!_     takes   active $
      w6j!_to!_w!_     takes inactive $
      wrk0 takes mode!_id!_ $
      end $
  if mode!_id!_ is Asitis then do
    begin $
      c!_to!_w3j!_     takes inactive $
      w3j!_to!_c!_     takes inactive $
      w!_to!_w6j!_     takes inactive $
      w6j!_to!_w!_     takes inactive $
      wrk0 takes mode!_id!_ $
      end $
  if mode!_id!_ is Show then do
    begin $
      write "Currently C.G. -> 3J is ", c!_to!_w3j!_, "." $
      write "Currently 3J -> C.G. is ", w3j!_to!_c!_, "." $
      write "Currently Racah -> 6J is ", w!_to!_w6j!_, "." $
      write "Currently 6J -> Racah is ", w6j!_to!_w!_, "." $
      end $
  if mode!_id!_ isnot Rose
     and
     mode!_id!_ isnot Wigner
     and
     mode!_id!_ isnot Asitis
                          then do
    begin $
      write "Usage: notation [ Rose, Wigner, Asitis, Show ]" $
      end $
  return wrk0 $
  end $
%
%
% If mode!_id!_ is
% abstract: the phases of the coefficients are left abstract,
% concrete: the phases of the coefficients are given in a
%           conventional manner. 
%
%              n * p!_                                 n
% Note: ( - 1 )         is the abstract form of ( - 1 )  .
%
procedure phasemode ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is abstract then phase!_ takes abstract $
    if mode!_id!_ is concrete then phase!_ takes concrete $
    if mode!_id!_ is show then write "Currently ", phase!_, "." $
    if mode!_id!_ isnot abstract and mode!_id!_ isnot concrete then
      write "Usage: phasemode [ abstract, concrete, show ]" $
    return phase!_ $
    end $
%
%
% If mode!_id!_ is
%   activate: factorials: factrl ( x + n ),
%             where x is a variable and n is an integer, 
%             are decomposed into the products of
%             a factorial:          factrl ( x )    and 
%             a partial factorial: pfactrl ( x, n ),
% deactivate: avoid factorial decomposition.  
%
procedure factred ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is deactivate then 
                                ftl!_red!_mode!_02 takes inactive $
    if mode!_id!_ is   activate then 
                                ftl!_red!_mode!_02 takes   active $
    if mode!_id!_ is       show then 
                      write "Currently ", ftl!_red!_mode!_02, "." $
    if mode!_id!_ isnot deactivate and mode!_id!_ isnot activate then
            write "Usage: factred [ deactivate, activate, show ]" $
    return ftl!_red!_mode!_02 $
    end $
%
%
% If mode!_id!_ is
%   activate: the ftcdiv!_ function is set ready to work,
%             where it directly transforms the ratios of 
%             two factorials: factrl ( x ) / factrl ( y )
%             into partial factorials: pfactrl ( y, x - y ),
% deactivate: avoid the ftcdiv!_ function to work. 
%
procedure factdiv ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is deactivate then 
                                ftl!_red!_mode!_03 takes inactive $
    if mode!_id!_ is   activate then 
                                ftl!_red!_mode!_03 takes   active $
    if mode!_id!_ is       show then 
                      write "Currently ", ftl!_red!_mode!_03, "." $
    if mode!_id!_ isnot deactivate and mode!_id!_ isnot activate then
            write "Usage: factred [ deactivate, activate, show ]" $
    return ftl!_red!_mode!_03 $
    end $
%
%
% If mode!_id!_ is
%   activate: the partial factorial functions: pfactrl ( x, n )
%             are replaced by actual pruduct of x + i 
%             ( i = 1, 2, ... n ):
%             ( x + 1 )( x + 2 )...( x + n ),
% deactivate: avoid this replacement. 
%              
procedure pfactred ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is   activate then 
                         reduction!_of!_pfactrl!_ takes   active $
    if mode!_id!_ is deactivate then 
                         reduction!_of!_pfactrl!_ takes inactive $
    if mode!_id!_ is       show then 
               write "Currently ", reduction!_of!_pfactrl!_, "." $
    if mode!_id!_ isnot activate and mode!_id!_ isnot deactivate then
      write "Usage: partialfactorial [ activate, deactivate, show ]" $
    ftl!_red!_mode!_01 takes reduction!_of!_pfactrl!_ $
    return reduction!_of!_pfactrl!_ $
    end $
%
%
% If mode!_id!_ is
%   activate: when both x and n are numbers, 
%             the partial factorial functions: pfactrl ( x, n )
%             are replaced by actual pruduct of x + i 
%             ( i = 1, 2, ... n ):
%             ( x + 1 )( x + 2 )...( x + n ),
% deactivate: avoid this replacement.
%             
% Note: The number 0 is replaced by a variable zero to avoid
%       a divide exception interruption by the system.
%       The expression "zero / zero" is replaced by 1 by the
%       system unconditionally. In this sense, the correctness
%       of the calculation is not guaranteed in general
%       when the present feature is on.
% 
procedure numericalpfactorial ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is reduce then 
                nonumerical!_pfactorial!_expression!_ takes   active $
    if mode!_id!_ is retain then 
                nonumerical!_pfactorial!_expression!_ takes inactive $
    if mode!_id!_ is   show then 
      write "Currently ", nonumerical!_pfactorial!_expression!_, "." $
    if mode!_id!_ isnot reduce and mode!_id!_ isnot retain then
         write "Usage: numericalpfactorial [ reduce, retain, show ]" $
    return nonumerical!_pfactorial!_expression!_ $
    end $
%
%
% If mode!_id!_ is
% reduce: replace inverse-of-factorial functions: tactrl ( x )
%         by inverse of factorial-functions:  1 / factrl ( x ),
% retain: leave tactrl ( x ) unchanged.
%
procedure inverseoffactorial ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is reduce then 
               noinverse!_of!_factorial!_expression!_ takes   active $
    if mode!_id!_ is retain then 
               noinverse!_of!_factorial!_expression!_ takes inactive $
    if mode!_id!_ is   show then 
      write "Currently ", noinverse!_of!_factorial!_expression!_, "."$
    if mode!_id!_ isnot reduce and mode!_id!_ isnot retain then
          write "Usage: inverseoffactorial [ reduce, retain, show ]" $
    return noinverse!_of!_factorial!_expression!_ $
    end $
%
%
% If mode!_id!_ is
% reduce: replace egg ( x ) by x for all x, 
% retain: avoid the replacement.
%
procedure chicksintheeggs ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is reduce then
               hatch!_egg!_mode!_ takes   active $
    if mode!_id!_ is retain then
               hatch!_egg!_mode!_ takes inactive $
    if mode!_id!_ is   show then
      write "Currently ", hatch!_egg!_mode!_, "."$
    if mode!_id!_ isnot reduce and mode!_id!_ isnot retain then
          write "Usage: chicksintheeggs [ reduce, retain, show ]" $
    return hatch!_egg!_mode!_ $
    end $
%
%
% If mode!_id!_ is
%   activate: set the following feature ready to work,       
% deactivate: avoid the following feature to work. 
%
% The feature is to replace the phase functions
% that came from the expression Gamma(x)Gamma(-x)
% by the ratios of polynomials. 
%
procedure polynomialformforphfunc ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is   activate then 
                     phase!_function!_red!_mode!_00!_ takes   active $
    if mode!_id!_ is deactivate then 
                     phase!_function!_red!_mode!_00!_ takes inactive $
    if mode!_id!_ is       show then 
           write "Currently ", phase!_function!_red!_mode!_00!_, "." $
    if mode!_id!_ isnot  activate and mode!_id!_ isnot deactivate then
      write 
     "Usage: polynomialformforphfunc [ activate, deactivate, show ]" $
    return phase!_function!_red!_mode!_00!_ $
    end $
%
%
% If mode!_id!_ is
%   activate: set the following feature ready to work,        
% deactivate: avoid the following feature to work.                      
%
% The feature is to replace the phase functions
% that came from the expression Gamma(x)Gamma(-x)
% by ( pi * x ) / sin ( pi * x ).          
%
% 
procedure sineformforphfunc ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is   activate then 
                     phase!_function!_red!_mode!_01!_ takes   active $
    if mode!_id!_ is deactivate then 
                     phase!_function!_red!_mode!_01!_ takes inactive $
    if mode!_id!_ is       show then 
           write "Currently ", phase!_function!_red!_mode!_01!_, "." $
    if mode!_id!_ isnot  activate and mode!_id!_ isnot deactivate then
      write 
           "Usage: sineformforphfunc [ activate, deactivate, show ]" $
    return phase!_function!_red!_mode!_01!_ $
    end $
%
%
% If mode!_id!_ is
%   activate: the system uses pfactrl ( x, n ) functions
%             as the intermediate form of the formula deduction. 
% deactivate: the system does not use pfactrl ( x, n ) functions   
%             as the intermediate form of the formula deduction.
%
% Note: This switch should normally be off. When the formula
%       reduction turned out to be very time consuming,
%       the user is encouraged to set this switch on
%       and repeat the reduction. He will be successful
%       if the God smiles at him.
%       
procedure indirectreduction ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is   activate then 
                         indirect!_reduction!_path!_ takes   active $
    if mode!_id!_ is deactivate then 
                         indirect!_reduction!_path!_ takes inactive $
    if mode!_id!_ is       show then 
               write "Currently ", indirect!_reduction!_path!_, "." $
    if mode!_id!_ isnot activate and mode!_id!_ isnot deactivate then
      write 
          "Usage: indirectreduction [ activate, deactivate, show ]" $
    return indirect!_reduction!_path!_ $
    end $
%
%
% If mode!_id!_ is
%   activate: the automatic internal reduction controlling feature
%             is set on,
% deactivate: the automatic internal reduction controlling feature
%             is set off.
%
% Note: The user is strongly recommended to set this switch on.
%       Otherwise, every manipulation must be done at his own risk.
%
procedure autocontrol ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is   activate then 
                         internal!_switch!_control!_ takes   active $
    if mode!_id!_ is deactivate then 
                         internal!_switch!_control!_ takes inactive $
    if mode!_id!_ is show then 
               write "Currently ", internal!_switch!_control!_, "." $
    if mode!_id!_ isnot activate and mode!_id!_ isnot deactivate then
          write "Usage: autocontrol [ activate, deactivate, show ]" $
    return internal!_switch!_control!_ $
    end $
%
%
% If mode!_id!_ is
%        final: specify to output the result in a conventional form,
% intermediate: specify to output the result in an intermediate form
%               that contains symbols internally used by the system.
%
% Note: The intermediate option is used to get the results 
%       that are usable in further reduce operations.
%
procedure outputform ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is final then do
      begin $
        numericalpfactorial reduce $
        inverseoffactorial reduce $
        chicksintheeggs reduce $
        out!_put!_form!_flag!_ takes final $
        end $
    if mode!_id!_ is intermediate then do
      begin $
        numericalpfactorial retain $
        inverseoffactorial retain $
        chicksintheeggs retain $
        out!_put!_form!_flag!_ takes intermediate $
        end $
    return out!_put!_form!_flag!_ $
    end $
%
%
% If mode!_id!_ is 
% chatty: run time interval reporting switch is set on,
%  quiet: run time interval reporting switch is set off.
%
procedure runtimereporting ( mode!_id!_ ) $
  begin $
    if mode!_id!_ is chatty then message!_level!_ takes chatty $
    if mode!_id!_ is quiet  then message!_level!_ takes quiet  $
    if mode!_id!_ is show then 
                     write "Currently ", message!_level!_, "." $
    if mode!_id!_ isnot chatty and mode!_id!_ isnot quiet then
       write "Usage: runtimereporting [ chatty, quiet, show ]" $
    return  message!_level!_ $
    end $
%
%
% Routines for arranging the runtime switches.
%
% The system has 4 runtime switches numbered 0, 1, 2, and 3.
%
%  Switch sw0!_ controls the factorial manipulation.
%  Switch sw1!_ controls the Clebsch-Gordan coefficient manipulation.
%  Switch sw2!_ controls the Racah coefficient manipulation.
%  Switch sw3!_ controls the X-coefficient manipulation.
%
%  Each switch has 8 dip switches. When a number 1 is specified at 
%  call of the following routines, the corresponding switch will be
%  set "active".
%  When a number 0 is specified at call of the following routines,
%  the corresponding switch will be set "inactive".
%  When any numbers or variables other than 1 or 0 is specified at 
%  call, the corresponding switch suffers no change.
%
%
% Set switch sw0!_.
%
procedure ftlctlsw 
  ( sw0!_0, sw0!_1, sw0!_2, sw0!_3, sw0!_4, sw0!_5, sw0!_6, sw0!_7 ) $
  begin scalar flag $
    flag takes unsuccessful $
    if sw0!_0 is 1 then ftl!_red!_mode!_00 takes   active $
    if sw0!_0 is 0 then ftl!_red!_mode!_00 takes inactive $
    if sw0!_1 is 1 then ftl!_red!_mode!_01 takes   active $
    if sw0!_1 is 0 then ftl!_red!_mode!_01 takes inactive $
    if sw0!_2 is 1 then ftl!_red!_mode!_02 takes   active $
    if sw0!_2 is 0 then ftl!_red!_mode!_02 takes inactive $
    if sw0!_3 is 1 then ftl!_red!_mode!_03 takes   active $
    if sw0!_3 is 0 then ftl!_red!_mode!_03 takes inactive $
    if sw0!_4 is 1 then ftl!_red!_mode!_04 takes   active $
    if sw0!_4 is 0 then ftl!_red!_mode!_04 takes inactive $
    if sw0!_5 is 1 then ftl!_red!_mode!_05 takes   active $
    if sw0!_5 is 0 then ftl!_red!_mode!_05 takes inactive $
    if sw0!_6 is 1 then ftl!_red!_mode!_06 takes   active $
    if sw0!_6 is 0 then ftl!_red!_mode!_06 takes inactive $
    if sw0!_7 is 1 then ftl!_red!_mode!_07 takes   active $
    if sw0!_7 is 0 then ftl!_red!_mode!_07 takes inactive $
    flag takes   successful $
    return flag $
    end $
%
procedure sw0!_ 
  ( sw0!_0, sw0!_1, sw0!_2, sw0!_3, sw0!_4, sw0!_5, sw0!_6, sw0!_7 ) $
  begin scalar flag $
    flag takes ftlctlsw 
  ( sw0!_0, sw0!_1, sw0!_2, sw0!_3, sw0!_4, sw0!_5, sw0!_6, sw0!_7 ) $
    return $
    end $
%
%
% Set switch sw1!_.
%
procedure cgcctlsw 
  ( sw1!_0, sw1!_1, sw1!_2, sw1!_3, sw1!_4, sw1!_5, sw1!_6, sw1!_7 ) $
  begin scalar flag $
    flag takes unsuccessful $
    if sw1!_0 is 1 then cgc!_red!_mode!_00 takes   active $
    if sw1!_0 is 0 then cgc!_red!_mode!_00 takes inactive $
    if sw1!_1 is 1 then cgc!_red!_mode!_01 takes   active $
    if sw1!_1 is 0 then cgc!_red!_mode!_01 takes inactive $
    if sw1!_2 is 1 then cgc!_red!_mode!_02 takes   active $
    if sw1!_2 is 0 then cgc!_red!_mode!_02 takes inactive $
    if sw1!_3 is 1 then cgc!_red!_mode!_03 takes   active $
    if sw1!_3 is 0 then cgc!_red!_mode!_03 takes inactive $
    if sw1!_4 is 1 then cgc!_red!_mode!_04 takes   active $
    if sw1!_4 is 0 then cgc!_red!_mode!_04 takes inactive $
    if sw1!_5 is 1 then cgc!_red!_mode!_05 takes   active $
    if sw1!_5 is 0 then cgc!_red!_mode!_05 takes inactive $
    if sw1!_6 is 1 then cgc!_red!_mode!_06 takes   active $
    if sw1!_6 is 0 then cgc!_red!_mode!_06 takes inactive $
    if sw1!_7 is 1 then cgc!_red!_mode!_07 takes   active $
    if sw1!_7 is 0 then cgc!_red!_mode!_07 takes inactive $
    flag takes   successful $
    return flag $
    end $
%
procedure sw1!_ 
  ( sw1!_0, sw1!_1, sw1!_2, sw1!_3, sw1!_4, sw1!_5, sw1!_6, sw1!_7 ) $
  begin scalar flag $
    flag takes cgcctlsw 
  ( sw1!_0, sw1!_1, sw1!_2, sw1!_3, sw1!_4, sw1!_5, sw1!_6, sw1!_7 ) $
    return $
    end $
%
%
% Set switch sw2!_.
%
procedure rccctlsw 
  ( sw2!_0, sw2!_1, sw2!_2, sw2!_3, sw2!_4, sw2!_5, sw2!_6, sw2!_7 ) $
  begin scalar flag $
    flag takes unsuccessful $
    if sw2!_0 is 1 then rcc!_red!_mode!_00 takes   active $
    if sw2!_0 is 0 then rcc!_red!_mode!_00 takes inactive $
    if sw2!_1 is 1 then rcc!_red!_mode!_01 takes   active $
    if sw2!_1 is 0 then rcc!_red!_mode!_01 takes inactive $
    if sw2!_2 is 1 then rcc!_red!_mode!_02 takes   active $
    if sw2!_2 is 0 then rcc!_red!_mode!_02 takes inactive $
    if sw2!_3 is 1 then rcc!_red!_mode!_03 takes   active $
    if sw2!_3 is 0 then rcc!_red!_mode!_03 takes inactive $
    if sw2!_4 is 1 then rcc!_red!_mode!_04 takes   active $
    if sw2!_4 is 0 then rcc!_red!_mode!_04 takes inactive $
    if sw2!_5 is 1 then rcc!_red!_mode!_05 takes   active $
    if sw2!_5 is 0 then rcc!_red!_mode!_05 takes inactive $
    if sw2!_6 is 1 then rcc!_red!_mode!_06 takes   active $
    if sw2!_6 is 0 then rcc!_red!_mode!_06 takes inactive $
    if sw2!_7 is 1 then rcc!_red!_mode!_07 takes   active $
    if sw2!_7 is 0 then rcc!_red!_mode!_07 takes inactive $
    flag takes   successful $
    return flag $
    end $
%
procedure sw2!_ 
  ( sw2!_0, sw2!_1, sw2!_2, sw2!_3, sw2!_4, sw2!_5, sw2!_6, sw2!_7 ) $
  begin scalar flag $
    flag takes rccctlsw 
  ( sw2!_0, sw2!_1, sw2!_2, sw2!_3, sw2!_4, sw2!_5, sw2!_6, sw2!_7 ) $
    return $
    end $
%
%
% Set switch sw3!_.
%
procedure xctctlsw 
  ( sw3!_0, sw3!_1, sw3!_2, sw3!_3, sw3!_4, sw3!_5, sw3!_6, sw3!_7 ) $
  begin scalar flag $
    flag takes unsuccessful $
    if sw3!_0 is 1 then xct!_red!_mode!_00 takes   active $
    if sw3!_0 is 0 then xct!_red!_mode!_00 takes inactive $
    if sw3!_1 is 1 then xct!_red!_mode!_01 takes   active $
    if sw3!_1 is 0 then xct!_red!_mode!_01 takes inactive $
    if sw3!_2 is 1 then xct!_red!_mode!_02 takes   active $
    if sw3!_2 is 0 then xct!_red!_mode!_02 takes inactive $
    if sw3!_3 is 1 then xct!_red!_mode!_03 takes   active $
    if sw3!_3 is 0 then xct!_red!_mode!_03 takes inactive $
    if sw3!_4 is 1 then xct!_red!_mode!_04 takes   active $
    if sw3!_4 is 0 then xct!_red!_mode!_04 takes inactive $
    if sw3!_5 is 1 then xct!_red!_mode!_05 takes   active $
    if sw3!_5 is 0 then xct!_red!_mode!_05 takes inactive $
    if sw3!_6 is 1 then xct!_red!_mode!_06 takes   active $
    if sw3!_6 is 0 then xct!_red!_mode!_06 takes inactive $
    if sw3!_7 is 1 then xct!_red!_mode!_07 takes   active $
    if sw3!_7 is 0 then xct!_red!_mode!_07 takes inactive $
    flag takes   successful $
    return flag $
    end $
%
procedure sw3!_ 
  ( sw3!_0, sw3!_1, sw3!_2, sw3!_3, sw3!_4, sw3!_5, sw3!_6, sw3!_7 ) $
  begin scalar flag $
    flag takes xctctlsw 
  ( sw3!_0, sw3!_1, sw3!_2, sw3!_3, sw3!_4, sw3!_5, sw3!_6, sw3!_7 ) $
    return $
    end $
%
%
% Set an individual dip switch separately.
%
procedure setswitch ( switch!_id!_, intnl!_no!_, sw!_position!_ ) $
  begin scalar flag $
    clear nop $
    dp0!_ takes nop $ dp1!_ takes nop $ 
    dp2!_ takes nop $ dp3!_ takes nop $
    dp4!_ takes nop $ dp5!_ takes nop $ 
    dp6!_ takes nop $ dp7!_ takes nop $
    if intnl!_no!_ is 0 then do
      begin $
        if sw!_position!_ is 1 then dp0!_ takes 1 $
        if sw!_position!_ is 0 then dp0!_ takes 0 $
        end $
    if intnl!_no!_ is 1 then do
      begin $
        if sw!_position!_ is 1 then dp1!_ takes 1 $
        if sw!_position!_ is 0 then dp1!_ takes 0 $
        end $
    if intnl!_no!_ is 2 then do
      begin $
        if sw!_position!_ is 1 then dp2!_ takes 1 $
        if sw!_position!_ is 0 then dp2!_ takes 0 $
        end $
    if intnl!_no!_ is 3 then do
      begin $
        if sw!_position!_ is 1 then dp3!_ takes 1 $
        if sw!_position!_ is 0 then dp3!_ takes 0 $
        end $
    if intnl!_no!_ is 4 then do
      begin $
        if sw!_position!_ is 1 then dp4!_ takes 1 $
        if sw!_position!_ is 0 then dp4!_ takes 0 $
        end $
    if intnl!_no!_ is 5 then do
      begin $
        if sw!_position!_ is 1 then dp5!_ takes 1 $
        if sw!_position!_ is 0 then dp5!_ takes 0 $
        end $
    if intnl!_no!_ is 6 then do
      begin $
        if sw!_position!_ is 1 then dp6!_ takes 1 $
        if sw!_position!_ is 0 then dp6!_ takes 0 $
        end $
    if intnl!_no!_ is 7 then do
      begin $
        if sw!_position!_ is 1 then dp7!_ takes 1 $
        if sw!_position!_ is 0 then dp7!_ takes 0 $
        end $
    if switch!_id!_ is 0 then do
      begin $
        flag takes sw0!_ 
          ( dp0!_, dp1!_, dp2!_, dp3!_, dp4!_, dp5!_, dp6!_, dp7!_ ) $
        end $
    if switch!_id!_ is 1 then do
      begin $
        flag takes sw1!_ 
          ( dp0!_, dp1!_, dp2!_, dp3!_, dp4!_, dp5!_, dp6!_, dp7!_ ) $
        end $
    if switch!_id!_ is 2 then do
      begin $
        flag takes sw2!_ 
          ( dp0!_, dp1!_, dp2!_, dp3!_, dp4!_, dp5!_, dp6!_, dp7!_ ) $
        end $
    if switch!_id!_ is 3 then do
      begin $
        flag takes sw3!_ 
          ( dp0!_, dp1!_, dp2!_, dp3!_, dp4!_, dp5!_, dp6!_, dp7!_ ) $
        end $
    return $
    end $
%
%
% Set the system defaults of the runtime switches.
%
procedure defaultswitchposition ( switch!_id!_ ) $
  begin scalar flag $
    if switch!_id!_ is 0 then do
      begin $
        flag takes sw0!_ ( 0, 0, 1, 1, 0, 0, 0, 0 ) $
        end $
    if switch!_id!_ is 1 then do
      begin $
        flag takes sw1!_ ( 1, 0, 0, 0, 0, 0, 0, 0 ) $
        end $
    if switch!_id!_ is 2 then do
      begin $
        flag takes sw2!_ ( 1, 0, 0, 0, 0, 0, 0, 0 ) $
        end $
    if switch!_id!_ is 3 then do
      begin $
        flag takes sw3!_ ( 1, 0, 0, 0, 0, 0, 0, 0 ) $
        end $
    return $
    end $
%
% Save current positions of the runtime switch sw0!_.
%
procedure save!_sw0!_ $
  begin $
    ftl!_red!_mode!_00!_save takes ftl!_red!_mode!_00 $
    ftl!_red!_mode!_01!_save takes ftl!_red!_mode!_01 $
    ftl!_red!_mode!_02!_save takes ftl!_red!_mode!_02 $
    ftl!_red!_mode!_03!_save takes ftl!_red!_mode!_03 $
    ftl!_red!_mode!_04!_save takes ftl!_red!_mode!_04 $
    ftl!_red!_mode!_05!_save takes ftl!_red!_mode!_05 $
    ftl!_red!_mode!_06!_save takes ftl!_red!_mode!_06 $
    ftl!_red!_mode!_07!_save takes ftl!_red!_mode!_07 $
    return $
    end $
%
% Restore the last positions of the runtime switch sw0!_.
%
% The call of this routine must be preceded by the call
% of the procedure: save!_sw0!_.
%
procedure restore!_sw0!_ $
  begin $
    ftl!_red!_mode!_00 takes ftl!_red!_mode!_00!_save $
    ftl!_red!_mode!_01 takes ftl!_red!_mode!_01!_save $
    ftl!_red!_mode!_02 takes ftl!_red!_mode!_02!_save $
    ftl!_red!_mode!_03 takes ftl!_red!_mode!_03!_save $
    ftl!_red!_mode!_04 takes ftl!_red!_mode!_04!_save $
    ftl!_red!_mode!_05 takes ftl!_red!_mode!_05!_save $
    ftl!_red!_mode!_06 takes ftl!_red!_mode!_06!_save $
    ftl!_red!_mode!_07 takes ftl!_red!_mode!_07!_save $
    return $
    end $
%
%
% Auxiliary functions to control the flow of the coefficient reduction.
% 
% The following features are effective when the flag: 
% internal!_switch!_control!_ is active.
%
procedure reduction!_control!_ ( cnt!_indx!_ ) $
  begin scalar flag $
    if internal!_switch!_control!_ is active and cnt!_indx!_ is 0
      then do
      begin $
        flag takes save!_sw0!_ () $
        end $
    if internal!_switch!_control!_ is active and cnt!_indx!_ is 1
      then do
      begin $
        if indirect!_reduction!_path!_ is active then do
          begin $
            flag takes sw0!_ 
            ( nop,   0, nop, nop, nop, nop, nop, nop ) $
            end
                                                 else do
          begin $
            if reduction!_of!_pfactrl!_    is active then do
              begin $
                flag takes sw0!_ 
                ( nop,   1, nop, nop, nop, nop, nop, nop ) $
                end
                                                     else do
              begin $
                flag takes sw0!_ 
                ( nop,   0, nop, nop, nop, nop, nop, nop ) $
                end $
            end $
        end $
    if internal!_switch!_control!_ is active and cnt!_indx!_ is 2
      then do
      begin $
        if reduction!_of!_pfactrl!_ is active then do
          begin $
            flag takes sw0!_ 
            ( nop,   1, nop, nop, nop, nop, nop, nop ) $
            end
                                              else do
          begin $
            flag takes sw0!_ 
            ( nop,   0, nop, nop, nop, nop, nop, nop ) $
            end $
        end $
    if internal!_switch!_control!_ is active  and cnt!_indx!_ is 3 
      then do
      begin $
        if noinverse!_of!_factorial!_expression!_ is active then do
          begin $
            flag takes sw0!_ 
            ( nop, nop, nop, nop, nop, nop, nop,   1 ) $
            end
                                                            else do
          begin $
            flag takes sw0!_ 
            ( nop, nop, nop, nop, nop, nop, nop,   0 ) $
            end $
        end $
    if internal!_switch!_control!_ is active  and cnt!_indx!_ is 4 
      then do
      begin $
        if nonumerical!_pfactorial!_expression!_ is active then do
          begin $
            flag takes sw0!_ 
            (   1, nop, nop, nop, nop, nop, nop, nop ) $
            end
                                                            else do
          begin $
            flag takes sw0!_ 
            (   0, nop, nop, nop, nop, nop, nop, nop ) $
            end $
        end $
    if internal!_switch!_control!_ is active  and cnt!_indx!_ is 5
      then do
      begin $
        if phase!_function!_red!_mode!_00!_ is active then do
          begin $
            flag takes sw0!_ 
            ( nop, nop, nop, nop,   1,   1, nop, nop ) $
            end
                                                            else do
          begin $
            flag takes sw0!_ 
            ( nop, nop, nop, nop, nop,   0, nop, nop ) $
            end $
        end $
    if internal!_switch!_control!_ is active  and cnt!_indx!_ is 6
      then do
      begin $
        if phase!_function!_red!_mode!_01!_ is active then do
          begin $
            flag takes sw0!_ 
            ( nop, nop, nop, nop,   1, nop,   1, nop ) $
            end
                                                            else do
          begin $
            flag takes sw0!_ 
            ( nop, nop, nop, nop, nop, nop,   0, nop ) $
            end $
        end $
    if internal!_switch!_control!_ is active and cnt!_indx!_ is 8
      then do
      begin $
        flag takes restore!_sw0!_ () $
        end $
    return flag $
    end $
%
%
%
% End of Prologue.  Define basic operators and mode setters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% Chapter 1. Collection of auxiliary operators and procedures          %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% The full path of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter10.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter10.r,v 1.1 92/04/16 20:37:40 koikef Locked $'/
%
%
%
% 1.1. factorial of a number or variable n:   n!
%      and its inverse:                     1/n!.
%
%  Notes (1) A phase variable p!_ must be either cleared out or set
%            equal to - 1 on purpose at any call of the following
%            functions that includes p!_. If it is set otherwise,
%            the result may be mathematically incorrect.
%        (2) The use of the following functions may be controlled at
%            reduction of the angular momentum coefficients.
%            To activate the usage, let the program switches
%            ftl!_red!_mode!_0# (# = 0, 1, 2, 3, 4, 5, 6, 7) be active.
%            To deactivate the usage let them be inactive.
%            The switch can be set by a call of ftlctlsw described
%            in Chapter 1.
%
%
operator  factrl   $
operator  factrl!_ $
operator  tactrl   $
operator  tactrl!_ $
operator pfactrl   $
operator pfactrl!_ $
operator ptactrl   $
operator ptactrl!_ $
operator egg       $
%
%
%  n! :
%
%
% Factorial of positive integer numbers n.
%
for all n such that fixp n and n > 0
  let factrl ( n ) be
  begin scalar k $
    factrl takes for k  takes 1 : n product k $
    return factrl $
    end $
for all n such that fixp n and n is 0
  let factrl ( n ) be 1 $
%
%
% Partial factorials for numerical arguments x and n.
%
% pfactrl ( x, n ) = (x+1)(x+2)(x+3).....(x+n).
%
for all x, n such that     fixp x and fixp n
                       and x >= 0 and n >=0
                       and ftl!_red!_mode!_00 is active
  let pfactrl ( x, n ) be factrl!_ ( x + n ) * tactrl!_ ( x ) $
for all x, n such that     fixp x and fixp n
                       and x <  0 and n >=0 and x + n + 1 <= 0
                       and ftl!_red!_mode!_00 is active
  let pfactrl ( x, n ) be 
    ( - 1 ) ** ( n * p!_ ) * pfactrl ( - x - n - 1, n ) $
for all x, n such that     fixp x and fixp n
                       and x <  0 and n >=0 and x + n + 1 >  0
                       and ftl!_red!_mode!_00 is active
  let pfactrl ( x, n ) be 
    ( - 1 ) ** ( ( - x - 1 ) * p!_ ) * pfactrl ( 0, - x - 1 )
    * zero * pfactrl ( 0, x + n ) $
for all x, n such that     fixp x and fixp n
                       and x >  0 and n < 0 and x + n >= 0
                       and ftl!_red!_mode!_00 is active
  let pfactrl ( x, n ) be ptactrl ( x + n, -n ) $
for all x, n such that     fixp x and fixp n
                       and x >  0 and n < 0 and x + n <  0
                       and ftl!_red!_mode!_00 is active
  let pfactrl ( x, n ) be 
    ( - 1 ) ** ( ( x + n + 1 ) * p!_ ) * ptactrl ( 0, - x - n - 1 )
    / zero * ptactrl ( 0, x ) $
for all x, n such that     fixp x and fixp n
                       and x <  0 and n < 0
                       and ftl!_red!_mode!_00 is active
  let pfactrl ( x, n ) be 
    ( - 1 ) ** ( n * p!_ ) * ptactrl ( - x - 1, - n ) $
%
%
% Partial factorials for a non-numerical first argument x.
%
% pfactrl ( x, n ) = (x+1)(x+2)(x+3).....(x+n).
%
for all x, n such that     not numberp x
                       and fixp n and n > 0
                       and ftl!_red!_mode!_01 is active
  let pfactrl ( x, n ) be
  begin scalar k, wrk1 $
    wrk1 takes for k takes 1 : n     product     egg ( x + k ) $
    return wrk1 $
    end $
%
for all x, n such that     not numberp x
                       and fixp n and n < 0
                       and ftl!_red!_mode!_01 is active
  let pfactrl ( x, n ) be
  begin scalar k, wrk1 $
    off exp $
    wrk1 takes for k takes n + 1 : 0 product 1 / egg ( x + k ) $
    return wrk1 $
    end $
%
for all x, y, n, m such that     not numberp x and not numberp y
                             and fixp n and fixp m
                             and y - x is n
  let pfactrl ( x, n ) * pfactrl ( y, m ) be pfactrl ( x, n + m ) $
%
for all x, y, n, m such that     not numberp x and not numberp y
                             and fixp n and fixp m
                             and x + y is 0
  let pfactrl ( x, n ) * pfactrl ( y, m ) be
    begin scalar wrk1 $
      if part ( x * dummy!_, 0 ) is minus then do
        begin $
          wrk1 takes   pfactrl ( y - n - 1, m + n + 1 )
                     / egg ( y ) * ( - 1 ) ** ( - n * p!_ ) $
          end
                                          else do
        begin $
          wrk1 takes   pfactrl ( x - m - 1, n + m + 1 )
                     / egg ( x ) * ( - 1 ) ** ( - m * p!_ ) $
          end $
      return wrk1 $
      end $
%
for all x, n such that fixp n and n is 0
  let pfactrl ( x, n ) be 1 $
%
for all x, n such that fixp n and n is 1
  let pfactrl ( x, n ) be egg ( x + n ) $
%
% Decompose factorial functions.
%
%  (x+n)! = x!(x+1)(x+2)...(x+n),
%  when  x is a variable and n is an integer number :
%
for all x, n such that     not numberp x
                       and fixp n and n > 0
                       and ftl!_red!_mode!_02 is active
  let factrl ( x + n ) be pfactrl ( x      ,   n ) * factrl!_ ( x ) $
%
for all x, n such that     not numberp x
                       and fixp n and n < 0
                       and ftl!_red!_mode!_02 is active
  let factrl ( x + n ) be ptactrl ( x + n  , - n ) * factrl!_ ( x ) $
%
for all x, n such that     not numberp x
                       and not fixp n 
                       and fixp ( 2 * n ) and 2 * n >   0
                       and ftl!_red!_mode!_02 is active
  let factrl ( x + n ) be 
    pfactrl ( x - 1/2,   n + 1/2 ) * factrl!_ ( x - 1/2 ) $
%
for all x, n such that     not numberp x
                       and not fixp n 
                       and fixp ( 2 * n ) and 2 * n < - 1
                       and ftl!_red!_mode!_02 is active
  let factrl ( x + n ) be 
    ptactrl ( x + n  , - n - 1/2 ) * factrl!_ ( x - 1/2 ) $
%
%
% Inverse of factorial of positive integer n.
%
%  1/n! :
%
for all n such that fixp n and n >= 0
  let tactrl ( n ) be
  begin scalar k $
    tactrl takes for k takes 1 : n product 1 / k $
    return tactrl $
    end $
%
for all n such that fixp n and n is 0
  let tactrl ( n ) be 1 $
%
%
%  1/(x+n)! = 1/{x!(x+1)(x+2)...(x+n)},
%  when  x is a variable and n is an integer number :
%
%
% Inverse of partial factorials for numerical arguments.
%
% ptactrl ( x, n ) = 1/{(x+1)(x+2)(x+3).....(x+n)}.
%
for all x, n such that     fixp x and fixp n
                       and x >= 0 and n >=0
                       and ftl!_red!_mode!_00 is active
  let ptactrl ( x, n ) be tactrl!_ ( x + n ) * factrl!_ ( x ) $
for all x, n such that     fixp x and fixp n
                       and x <  0 and n >=0 and x + n + 1 <= 0
                       and ftl!_red!_mode!_00 is active
  let ptactrl ( x, n ) be 
    ( - 1 ) ** ( - n * p!_ ) * ptactrl ( - x - n - 1, n ) $
for all x, n such that     fixp x and fixp n
                       and x <  0 and n >=0 and x + n + 1 >  0
                       and ftl!_red!_mode!_00 is active
  let ptactrl ( x, n ) be 
    ( - 1 ) ** ( ( x + 1 ) * p!_ ) * ptactrl ( 0, - x - 1 )
    / zero * ptactrl ( 0, x + n ) $
for all x, n such that     fixp x and fixp n
                       and x >  0 and n < 0 and x + n >= 0
                       and ftl!_red!_mode!_00 is active
  let ptactrl ( x, n ) be pfactrl ( x + n, -n ) $
for all x, n such that     fixp x and fixp n
                       and x >  0 and n < 0 and x + n <  0
                       and ftl!_red!_mode!_00 is active
  let ptactrl ( x, n ) be 
    ( - 1 ) ** ( - ( x + n + 1 ) * p!_ ) * pfactrl ( 0, - x - n - 1 )
    * zero * pfactrl ( 0, x ) $
for all x, n such that     fixp x and fixp n
                       and x <  0 and n < 0
                       and ftl!_red!_mode!_00 is active
  let ptactrl ( x, n ) be 
    ( - 1 ) ** ( - n * p!_ ) * pfactrl ( - x - 1, - n ) $
%
%
% Inverse of partial factorials for a non-numerical first argument.
%
% ptactrl ( x, n ) = 1/{(x+1)(x+2)(x+3).....(x+n)}.
%
for all x, n such that     not numberp x and fixp n and n > 0
                       and ftl!_red!_mode!_01 is active
  let ptactrl ( x, n ) be
  begin scalar k, wrk1 $
    off exp $
    wrk1 takes for k takes 1 : n     product 1 / egg ( x + k ) $
    return wrk1 $
    end $
%
for all x, n such that     not numberp x and fixp n and n < 0
                       and ftl!_red!_mode!_01 is active
  let ptactrl ( x, n ) be
  begin scalar k, wrk1 $
    off exp $
    wrk1 takes for k takes n + 1 : 0 product     egg ( x + k ) $
    return wrk1 $
    end $
%
for all x, y, n, m such that     not numberp x and not numberp y
                             and fixp n and fixp m
                             and y - x is n
  let ptactrl ( x, n ) * ptactrl ( y, m ) be ptactrl ( x, n + m ) $
%
for all x, y, n, m such that     not numberp x and not numberp y
                             and fixp n and fixp m
                             and x + y is 0
  let ptactrl ( x, n ) * ptactrl ( y, m ) be
    begin scalar wrk1 $
      if part ( x * dummy!_, 0 ) is minus then do
        begin $
          wrk1 takes   ptactrl ( y - n - 1, m + n + 1 )
                     * egg ( y ) / ( - 1 ) ** ( - n * p!_ ) $
          end
                                          else do
        begin $
          wrk1 takes   ptactrl ( x - m - 1, n + m + 1 )
                     * egg ( x ) / ( - 1 ) ** ( - m * p!_ ) $
          end $
      return wrk1 $
      end $
%
for all x, n such that fixp n and n is 0
  let ptactrl ( x, n ) be 1 $
%
for all x, n such that fixp n and n is 1
  let ptactrl ( x, n ) be 1 / egg ( x + n ) $
%
%
% Decompose inverse of factorial functions.
%
%  1/(x+n)! = 1/{x!(x+1)(x+2)...(x+n)},
%  when  x is a variable and n is an integer number :
%
for all x, n such that     not numberp x
                       and fixp n and n > 0
                       and ftl!_red!_mode!_02 is active
  let tactrl ( x + n ) be ptactrl ( x      ,   n ) * tactrl!_ ( x ) $
%
for all x, n such that     not numberp x
                       and fixp n and n < 0
                       and ftl!_red!_mode!_02 is active
  let tactrl ( x + n ) be pfactrl ( x + n  , - n ) * tactrl!_ ( x ) $
%
for all x, n such that     not numberp x
                       and not fixp n 
                       and fixp ( 2 * n ) and 2 * n >   0
                       and ftl!_red!_mode!_02 is active
  let tactrl ( x + n ) be 
    ptactrl ( x - 1/2,   n + 1/2 ) * tactrl!_ ( x - 1/2 ) $
%
for all x, n such that     not numberp x
                       and not fixp n 
                       and fixp ( 2 * n ) and ( 2 * n ) < - 1
                       and ftl!_red!_mode!_02 is active
  let tactrl ( x + n ) be 
    pfactrl ( x + n  , - n - 1/2 ) * tactrl!_ ( x - 1/2 ) $
%
%
% Contract the ratio of factorial functions.
%
%  x!/y! = x(x-1)...(y+2)(y+1),
%  when x-y is an integer number.
%
operator ftcdiv!_ $
%
for all x, y such that     not fixp ( x - y )
                       or  ( numberp x and numberp y )
                       or  ftl!_red!_mode!_03 is inactive
  let ftcdiv!_ ( x, y ) be factrl!_ ( x ) * tactrl!_ ( y ) $
%
for all x, y such that     fixp ( x - y ) and x - y isnot 0
                       and ( not numberp x or not numberp y )
                       and ftl!_red!_mode!_03 is active
  let ftcdiv!_ ( x, y ) be
    begin scalar wrk1 $
      if x - y > 0 then wrk1 takes pfactrl ( y, x - y ) $
      if x - y < 0 then wrk1 takes ptactrl ( x, y - x ) $
      return wrk1 $
      end $
%
for all x, y such that     fixp ( x - y ) and x - y is 0
  let ftcdiv!_ ( x, y ) be 1 $
%
for all x, y such that     fixp ( x - y ) and x - y isnot 0
                       and ( not numberp x or not numberp y )
                       and ftl!_red!_mode!_03 is active
  let factrl ( x ) * tactrl ( y ) be
    begin scalar wrk1 $
      if x - y > 0 then wrk1 takes pfactrl ( y, x - y ) $
      if x - y < 0 then wrk1 takes ptactrl ( x, y - x ) $
      return wrk1 $
      end $
%
for all x, y such that     fixp ( x - y ) and x - y is 0
  let factrl ( x ) * tactrl ( y ) be 1 $
%
for all x, y, n, m such that     fixp ( x - y ) and fixp n and fixp m
                             and x - y < m and x - y > 0
  let pfactrl ( x, n ) * ptactrl ( y, m ) be
    begin scalar wrk1 $
      if x - y > m - n then wrk1 takes
         ptactrl ( y, x - y ) * pfactrl ( y + m, x - y + n - m ) $
      if x - y = m - n then wrk1 takes
         ptactrl ( y, x - y ) $
      if x - y < m - n then wrk1 takes
         ptactrl ( y, x - y ) * ptactrl ( x + n, y - x + m - n ) $
      return wrk1 $
      end $
%
for all x, y, n, m such that     fixp ( y - x ) and fixp n and fixp m
                             and y - x < n and y - x > 0
  let pfactrl ( x, n ) * ptactrl ( y, m ) be
    begin scalar wrk1 $
      if y - x > n - m then wrk1 takes
         pfactrl ( x, y - x ) * ptactrl ( x + n, y - x + m - n ) $
      if y - x = n - m then wrk1 takes
         pfactrl ( x, y - x ) $
      if y - x < n - m then wrk1 takes
         pfactrl ( x, y - x ) * pfactrl ( y + m, x - y + n - m ) $
      return wrk1 $
      end $
%
for all x, y, n, m such that     fixp ( x - y ) and x - y is    0
                             and fixp ( n - m )
  let pfactrl ( x, n ) * ptactrl ( y, m ) be
    begin scalar wrk1 $
      if n - m > 0 then wrk1 takes pfactrl ( x + m, n - m ) $
      if n - m = 0 then wrk1 takes 1 $
      if n - m < 0 then wrk1 takes ptactrl ( x + n, m - n ) $
      return wrk1 $
      end $
%
%
% Contract the product of factorial functions.
%
%  x!y! = x!((x+y)-x)! = pi*x/sin(pi*x)(1-x)(2-x)...(y-1)y
%  when x+y is an integer number.
%
%  Note tat x!(-x)! = pi*x/sin(pi*x) in general.
%
operator phase!_function!_ $
%
for all x, y such that     x + y is 0 and not numberp x
                       and ftl!_red!_mode!_04 is active
  let factrl ( x ) * factrl ( y ) be
    begin scalar wrk1 $
  if part ( x * dummy!_, 0 ) is minus then do
    begin $
      wrk1 takes phase!_function!_ ( y ) $
      end
                                      else do
    begin $
      wrk1 takes phase!_function!_ ( x ) $
      end $
  return wrk1 $
  end $
%
for all x, y such that     x + y is - 1 and not numberp x
                       and ftl!_red!_mode!_04 is active
  let factrl ( x ) * factrl ( y ) be
    begin scalar wrk1 $
  if part ( x * dummy!_, 0 ) is minus then do
    begin $
      wrk1 takes phase!_function!_ ( y ) / egg ( x + 1 ) $
      end
                                      else do
    begin $
      wrk1 takes phase!_function!_ ( x ) / egg ( y + 1 ) $
      end $
  return wrk1 $
  end $
%
for all x such that     ftl!_red!_mode!_05 is active
                    and check!_var!_range!_ ( x ) is determined
  let phase!_function!_ ( x ) be
    begin scalar jidx, kidx, start, stop, wrk0, wrk1, wrk2 $
      start takes lower!_bound!_of!_x!_ $
      stop  takes upper!_bound!_of!_x!_ $
      kidx takes start $
      wrk1 takes 0 $
      while kidx <= stop do
        begin $
          jidx takes start $
          wrk0 takes 1 $
          while jidx <= stop do
            begin $
              on exp $
              if jidx isnot kidx then wrk0 takes wrk0 * ( x - jidx ) $
              jidx takes jidx + 1 $
              end $
          on exp $
          wrk1 takes wrk1 + ( - 1 ) ** ( - kidx ) * wrk0 $
          kidx takes kidx + 1 $
          end $
      wrk2 takes   ptactrl( x - stop - 1, stop - start + 1 ) 
                 * egg ( wrk1 ) * egg ( x ) $
      return wrk2 $
      end $
%
for all x such that ftl!_red!_mode!_06 is active
  let phase!_function!_ ( x ) be pi * egg ( x ) / sin ( pi * x ) $
%
%
% Decompose the inverse of the product of factorials.
%
%  1/x!y! = 1/x!((x+y)-x)! = 1/{pi*x/sin(pi*x)(1-x)(2-x)...(y-1)y}
%  when x+y is an integer number. 
%
%  Note that x!(-x)! = pi*x/sin(pi*x) in general.
%
operator qhase!_function!_ $
%
for all x, y such that     x + y is 0 and not numberp x
                       and ftl!_red!_mode!_04 is active
  let tactrl ( x ) * tactrl ( y ) be
    begin scalar wrk1 $
      if part ( x * dummy!_, 0 ) is minus then do
        begin $
          wrk1 takes qhase!_function!_ ( y ) $
          end
                                          else do
        begin $
          wrk1 takes qhase!_function!_ ( x ) $
          end $
      return wrk1 $
      end $
%
for all x, y such that     x + y is - 1 and not numberp x
                       and ftl!_red!_mode!_04 is active
  let tactrl ( x ) * tactrl ( y ) be
    begin scalar wrk1 $
      if part ( x * dummy!_, 0 ) is minus then do
        begin $
          wrk1 takes qhase!_function!_ ( y ) * egg ( x + 1 ) $
          end
                                          else do
        begin $
          wrk1 takes qhase!_function!_ ( x ) * egg ( y + 1 ) $
          end $
      return wrk1 $
      end $
%
for all x such that     ftl!_red!_mode!_05 is active
                    and check!_var!_range!_ ( x ) is determined
  let qhase!_function!_ ( x ) be
    begin scalar jidx, kidx, start, stop, wrk0, wrk1, wrk2 $
      start takes lower!_bound!_of!_x!_ $
      stop  takes upper!_bound!_of!_x!_ $
      kidx takes start $
      wrk1 takes 0 $
      while kidx <= stop do
        begin $
          jidx takes start $
          wrk0 takes 1 $
          while jidx <= stop do
            begin $
              on exp $
              if jidx isnot kidx then wrk0 takes wrk0 * ( x - jidx ) $
              jidx takes jidx + 1 $
              end $
          on exp $
          wrk1 takes wrk1 + ( - 1 ) ** ( - kidx  ) * wrk0 $
          kidx takes kidx + 1 $
          end $
      wrk2 takes   pfactrl( x - stop - 1, stop - start + 1 ) 
                 / egg ( wrk1 ) / egg ( x ) $
      return wrk2 $
      end $
%
for all x such that ftl!_red!_mode!_06 is active
  let qhase!_function!_ ( x ) be sin ( pi * x ) / ( pi * egg ( x ) ) $
%
%
for all x, y such that x - y is 0
  let phase!_function!_ ( x ) * qhase!_function!_ ( y ) be 1 $
%
%
for all x, y such that     x + y is 0 and not numberp x
                       and ftl!_red!_mode!_04 is active
  let factrl ( x ) * tactrl ( y ) be
    begin scalar wrk1 $
  if part ( x * dummy!_, 0 ) is minus then do
    begin $
      wrk1 takes phase!_function!_ ( y ) * tactrl!_ ( y ) ** 2 $
      end
                                      else do
    begin $
      wrk1 takes qhase!_function!_ ( x ) * factrl!_ ( x ) ** 2 $
      end $
  return wrk1 $
  end $
%
for all x, y such that     x + y is - 1 and not numberp x
                       and ftl!_red!_mode!_04 is active
  let factrl ( x ) * tactrl ( y ) be
    begin scalar wrk1 $
  if part ( x * dummy!_, 0 ) is minus then do
    begin $
      wrk1 takes phase!_function!_ ( y ) 
                 * ptactrl ( x, 1 ) * tactrl!_ ( y ) ** 2 $
      end
                                      else do
    begin $
      wrk1 takes qhase!_function!_ ( x ) 
                 * pfactrl ( y, 1 ) * factrl!_ ( x ) ** 2 $
      end $
  return wrk1 $
  end $
%
%
%  The following instruction is used to eliminate
%  the expression "tactrl" from the final output.
%
%
for all x such that ftl!_red!_mode!_07 is active
  let tactrl ( x ) be 1 / factrl ( x ) $
%
for all x, n  such that ftl!_red!_mode!_07 is active
  let ptactrl ( x, n ) be 1 / pfactrl ( x, n ) $
%
%
%  Gateway routines used for making smooth and quick
%  reductions of the factorial expressions.
%
%
for all x
  let factrl!_ ( x ) be
    begin scalar wrk0, wrk1, wrk2 $
      on exp $
      off mcd $
      wrk0 takes x $
      wrk1 takes factrl ( wrk0 ) $
      on mcd $
      wrk2 takes wrk1 $
      return wrk2 $
      end $
%
%
for all x
  let tactrl!_ ( x ) be
    begin scalar wrk0, wrk1, wrk2 $
      on exp $
      off mcd $
      wrk0 takes x $
      wrk1 takes tactrl ( wrk0 ) $
      on mcd $
      wrk2 takes wrk1 $
      return wrk2 $
      end $
%
%
% 1.2. Kroneker's delta:
%
%
operator krdlta   $
operator krdlta!_ $
%
for all x, y such that x  is  y
  let krdlta ( x, y ) be 1 $
%
for all x, y such that numberp ( x - y ) and ( x - y ) neq 0
  let krdlta ( x, y ) be 0 $
%
for all x, y
  let krdlta ( x, y ) * krdlta ( y, x ) be krdlta ( x, y ) $
%
for all x, y, n such that not numberp x
  let krdlta ( x + n, y ) be krdlta ( x, y - n ) $
%
for all x, y, n such that not numberp x
  let krdlta ( x * n, y ) be krdlta ( x, y / n ) $
%
for all x, y such that numberp x and not numberp y
  let krdlta ( x, y ) be krdlta ( y, x ) $
%
for all x, y
  let krdlta ( x, y ) ** 2 be krdlta ( x, y ) $
%
%
% The following routine is prepared on purpose of substituting
% y into x in the expressions containing krdlta ( x, y ).
%
for all x, y such that not numberp x and x isnot y
                       and meta!_mode!_of!_krdlta!_ is active
                       and lterm ( x, mainvar x ) is x
  let krdlta ( x, y ) be
    begin scalar wrk0 $
      meta!_mode!_of!_krdlta!_ takes inactive $
      wrk0 takes sub ( x = y, work!_for!_meta!_mode!_of!_krdlta!_ ) $
      work!_for!_meta!_mode!_of!_krdlta!_ takes 
                                           krdlta!_ ( x, y ) * wrk0 $
      meta!_mode!_of!_krdlta!_ takes active $
      return 1 $
      end $
%
for all x, y such that meta!_mode!_of!_krdlta!_ is inactive
  let krdlta!_ ( x, y ) be krdlta ( x, y ) $
%
%
% 1.3. Step function:
%
%
operator stepup $
%
for all x, y such that x is y
  let stepup ( x, y ) be 1 $
%
for all x, y such that fixp ( x - y ) and ( x - y ) >= 0
  let stepup ( x, y ) be 1 $
%
for all x, y such that fixp ( x - y ) and ( x - y ) <  0
  let stepup ( x, y ) be 0 $
%
for all x, y, z such that numberp ( y - z ) and y - z >= 0
  let stepup ( x, y ) * stepup ( x, z ) be stepup ( x, y ) $
%
for all x, y
  let stepup ( x, y ) ** 2 be stepup ( x, y ) $
%
%
operator stepdown $
%
for all x, y such that x is y
  let stepdown ( x, y ) be 1 $
%
for all x, y such that fixp ( x - y ) and ( x - y ) >  0
  let stepdown ( x, y ) be 0 $
%
for all x, y such that fixp ( x - y ) and ( x - y ) <= 0
  let stepdown ( x, y ) be 1 $
%
for all x, y, z such that numberp ( y - z ) and y - z <= 0
  let stepdown ( x, y ) * stepdown ( x, z ) be stepdown ( x, y ) $
%
for all x, y
  let stepdown ( x, y ) ** 2 be stepdown ( x, y ) $
%
for all x, y
  let sqrt ( stepdown ( x, y ) ) be stepdown ( x, y ) $
%
%
operator window   $
operator window!_ $
%
for all x, y, z such that x is y or x is z
  let window ( x, y, z ) be 1 $
%
for all x, y, z such that     numberp ( x - y ) and x - y >= 0
                          and numberp ( x - z ) and x - z <= 0
  let window ( x, y, z ) be 1 $
%
for all x, y, z such that     numberp ( x - y ) and x - y <  0
                           or numberp ( x - z ) and x - z >  0
  let window ( x, y, z ) be 0 $
%
for all x, y, z such that     numberp ( y - z ) and y - z >  0
  let window ( x, y, z ) be 0 $
%
for all x, y, z such that     numberp ( y - z ) and y - z =  0
  let window ( x, y, z ) be krdlta ( x, y ) $
%
for all x, y, z, n such that numberp n
  let window ( x + n, y, z ) be window!_ ( x, y - n, z - n ) $
%
for all x, y1, z1, y2, z2 
  such that     numberp ( y1 - y2 ) and y1 - y2 >= 0
            and numberp ( z1 - z2 ) and z1 - z2 <= 0
  let window ( x, y1, z1 ) * window ( x, y2, z2 ) be 
    window!_ ( x, y1, z1 ) $
%
for all x, y1, z1, y2, z2 
  such that     numberp ( y1 - y2 ) and y1 - y2 >= 0
            and numberp ( z1 - z2 ) and z1 - z2 >  0
  let window ( x, y1, z1 ) * window ( x, y2, z2 ) be 
    window!_ ( x, y1, z2 ) $
%
for all x, y1, z1, y2, z2 
  such that     numberp ( y1 - y2 ) and y1 - y2 <  0
            and numberp ( z1 - z2 ) and z1 - z2 >  0
  let window ( x, y1, z1 ) * window ( x, y2, z2 ) be 
    window!_ ( x, y2, z2 ) $
%
for all x, y1, z1, y2, z2 
  such that     numberp ( y1 - y2 ) and y1 - y2 <  0
            and numberp ( z1 - z2 ) and z1 - z2 <= 0
  let window ( x, y1, z1 ) * window ( x, y2, z2 ) be 
    window!_ ( x, y2, z1 ) $
%
for all x, y, z
  let window ( x, y, z ) ** 2 be window!_ ( x, y, z ) $
%
for all x, y, z
  let window!_ ( x, y, z ) be
    begin scalar wrk1, xwrk, ywrk, zwrk $
      on exp $ off mcd $
      xwrk takes x $ ywrk takes y $ zwrk takes z $
      if x isnot 0 and part ( dummy1!_ * dummy2!_ * xwrk, 0 ) is minus
        then wrk1 takes window ( - xwrk, - zwrk, - ywrk )
        else wrk1 takes window (   xwrk,   ywrk,   zwrk ) $
      off exp $ on mcd $
      return wrk1 $
      end $
%
%
% The following procedure creates a window function
% for a given argument according to the range of
% its variation.
%
procedure make!_window!_for!_x!_ ( x ) $
  begin scalar wrk1 $
    wrk1 takes x $
    if not numberp x then do
      begin $
        if check!_var!_range!_ ( x ) is determined then do
          begin $
            ww!_x!_ takes 
              ww!_x!_ * window!_ ( x, lower!_bound!_of!_x!_, 
                                      upper!_bound!_of!_x!_ ) $
            if lower!_bound!_of!_x!_ is upper!_bound!_of!_x!_ then
              wrk1 takes lower!_bound!_of!_x!_ $
            end $
        end $
    return wrk1 $
    end $
%
%
% The following set of procedures determines
% the range of angular momentum variable
% according to the given triangular conditions
% and the condition that the magnetic component
% of an angular momentum must not exceed its own length.
%
procedure check!_var!_range!_aux!_ ( x, y ) $
  begin scalar wrk1, wrk2, wrk3, wrk4, wrk5, 
        wrk6, wrk7, wrk8, wrk9, flag $
    flag takes undetermined $
    if not numberp x and not numberp y then do
      begin $
        if fixp ( 2 * ( y - x ) ) then do
          begin $
            wrk1 takes - y + x $
            if     lower!_bound!_of!_x!_ is undetermined
                or fixp lower!_bound!_of!_x!_
               and lower!_bound!_of!_x!_< wrk1
               then lower!_bound!_of!_x!_ takes wrk1 $
            if fixp ( y - x ) then do
              begin $
                if type!_of!_var!_x!_ is undetermined then 
                   type!_of!_var!_x!_ takes  fullinteger $
                if type!_of!_var!_x!_ is   halfinteger then 
                   type!_of!_var!_x!_ takes typeconflict $
                end
                              else do
              begin $
                if type!_of!_var!_x!_ is undetermined then 
                   type!_of!_var!_x!_ takes  halfinteger $
                if type!_of!_var!_x!_ is  fullinteger then
                   type!_of!_var!_x!_ takes typeconflict $
                end $
            end
                          else
        if fixp ( 2 * ( y + x ) ) then do
          begin $
            wrk1 takes   y + x $
            if     upper!_bound!_of!_x!_ is undetermined
                or fixp upper!_bound!_of!_x!_
               and upper!_bound!_of!_x!_> wrk1
               then upper!_bound!_of!_x!_ takes wrk1 $
            if fixp ( y + x ) then do
              begin $
                if type!_of!_var!_x!_ is undetermined then 
                   type!_of!_var!_x!_ takes  fullinteger $
                if type!_of!_var!_x!_ is  halfinteger then
                   type!_of!_var!_x!_ takes typeconflict $
                end
                              else do
              begin $
                if type!_of!_var!_x!_ is undetermined then 
                   type!_of!_var!_x!_ takes  halfinteger $
                if type!_of!_var!_x!_ is  fullinteger then
                   type!_of!_var!_x!_ takes typeconflict $
                end $
            end
                          else do
          begin $
            on exp $
            off mcd $
            wrk2 takes mainvar x $
            wrk3 takes lterm ( x, wrk2 ) $
            wrk4 takes lcof  ( x, wrk2 ) $
            wrk5 takes x - wrk3 $
            wrk6 takes mainvar y $
            wrk7 takes lterm ( y, wrk6 ) $
            wrk8 takes lcof  ( y, wrk6 ) $
            wrk9 takes y - wrk7 $
            off exp $
            on mcd $
            if     wrk2 is wrk6
               and wrk3 is wrk2 * wrk4 and numberp wrk5
               and wrk7 is wrk6 * wrk8 and numberp wrk9
               then do
                 begin $
                   wrk1 takes wrk4 * ( - wrk9 / wrk8 ) + wrk5 $
                   if wrk4 * wrk8 > 0 then do
                     begin $
                       if     lower!_bound!_of!_x!_ is undetermined
                           or fixp lower!_bound!_of!_x!_
                          and lower!_bound!_of!_x!_< wrk1
                          then lower!_bound!_of!_x!_ takes wrk1 $
                       end $
                   if wrk4 * wrk8 < 0 then do
                     begin $
                       if     upper!_bound!_of!_x!_ is undetermined
                           or fixp upper!_bound!_of!_x!_
                          and upper!_bound!_of!_x!_> wrk1
                          then upper!_bound!_of!_x!_ takes wrk1 $
                       end $
                   end $
            end $
        end $
    if     upper!_bound!_of!_x!_ isnot undetermined
       and lower!_bound!_of!_x!_ isnot undetermined
       then flag takes determined $
    return flag $
    end $
%
procedure check!_var!_range!_ ( x ) $
  begin scalar wrk1, wrk2, flag $
    if numberp x then do
      begin $
        lower!_bound!_of!_x!_ takes x $
        upper!_bound!_of!_x!_ takes x $
        if fixp ( 2 * x + 1 ) 
          then type!_of!_var!_x!_ takes  halfinteger
          else type!_of!_var!_x!_ takes typeconflict $
        if fixp       x       
          then type!_of!_var!_x!_ takes  fullinteger $
        wrk2 takes determined $
        end
              else do
      begin $
        lower!_bound!_of!_x!_ takes undetermined $
        upper!_bound!_of!_x!_ takes undetermined $
        type!_of!_var!_x!_ takes undetermined $
          if triangular!_cond!_check!_mode!_ is Xcoef then do
            begin $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_1!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_2!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_3!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_4!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_5!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_6!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_7!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_8!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_9!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_a!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_b!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_c!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_d!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_e!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_f!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_g!_save!_ ) $
              flag takes check!_var!_range!_aux!_ 
                ( x , Xcoef!_trang!_cond!_h!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_i!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_j!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_k!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_l!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_m!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_n!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_o!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_p!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_q!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Xcoef!_trang!_cond!_r!_save!_ ) $
              end $
          if triangular!_cond!_check!_mode!_ is Racah then do
            begin $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_1!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_2!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_3!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_4!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_5!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_6!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_7!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_8!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_9!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_a!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_b!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_c!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_d!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_e!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_f!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_g!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_h!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Racah!_trang!_cond!_i!_save!_ ) $
              end $
          if triangular!_cond!_check!_mode!_ is Clebsch then do
            begin $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_1!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_2!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_3!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_4!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_5!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_6!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_7!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_8!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_9!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_a!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_b!_save!_ ) $
              flag takes check!_var!_range!_aux!_
                ( x , Clebsch!_trang!_cond!_c!_save!_ ) $
              end $
        if     lower!_bound!_of!_x!_ is undetermined
            or upper!_bound!_of!_x!_ is undetermined
            or numberp 
                 ( lower!_bound!_of!_x!_ - upper!_bound!_of!_x!_ )
           and lower!_bound!_of!_x!_ - upper!_bound!_of!_x!_ > 0
           then wrk2 takes undetermined
           else wrk2 takes   determined $
        if wrk2 is determined and type!_of!_var!_x!_ is fullinteger 
          then do
          begin $
            if     not fixp lower!_bound!_of!_x!_
               and fixp ( 2 * lower!_bound!_of!_x!_ + 1 ) then
                   lower!_bound!_of!_x!_ takes 
                                         lower!_bound!_of!_x!_ + 1/2 $
            if     not fixp upper!_bound!_of!_x!_
               and fixp ( 2 * upper!_bound!_of!_x!_ + 1 ) then
                   upper!_bound!_of!_x!_ takes 
                                         upper!_bound!_of!_x!_ - 1/2 $
            end $
        if wrk2 is determined and type!_of!_var!_x!_ is halfinteger 
          then do
          begin $
            if fixp lower!_bound!_of!_x!_ then
              lower!_bound!_of!_x!_ takes 
                                    lower!_bound!_of!_x!_ + 1/2 $
            if fixp upper!_bound!_of!_x!_ then
              upper!_bound!_of!_x!_ takes 
                                    upper!_bound!_of!_x!_ - 1/2 $
            end $
        end $
    return wrk2 $
    end $
%
%
% A procedure to replace the "egg" function by its own argument.
%
% The mode controller hatch!_egg!_mode!_ can be set by calling:
% "chicksintheeggs ( reduce or retain )".
% Here the reduce option enables this replacement, and the retain
% option disables it.
%
procedure hatcheggs ( x ) $
    begin scalar wrk1 $
      if hatch!_egg!_mode!_ is active then do
        begin $
          off exp $
          for all y let egg ( y ) be y $
          wrk1 takes x $
          for all y clear egg ( y ) $
          off exp $
          end 
                                      else do
        begin $
          off exp $
          wrk1 takes x $
          off exp $
          end $ 
      return wrk1 $
    end $
%
%
% A procedure to replace the half egg function: sqrt ( egg x )  
% by a simple sqrt function: sqrt ( x ).
%
% Note: The system use only. 
%       The call of this procedure from the top level
%       could be errorness.
%
operator halfegg!_ $
%
procedure hatch!_half!_eggs!_ ( x ) $
  begin scalar wrk1, wrk2 $
    on exp $
    wrk1 takes x $
    for all y    let sqrt ( egg ( y ) ) 
      be halfegg!_ ( y )              $
    for all y, z let sqrt ( egg ( y ) * z )
      be halfegg!_ ( y ) * sqrt ( z ) $
    for all y, z let sqrt ( egg ( y ) ) * z 
      be halfegg!_ ( y ) *        z   $
    wrk2 takes wrk1 $
    for all y    clear sqrt ( egg ( y ) )     $
    for all y, z clear sqrt ( egg ( y ) * z ) $
    for all y, z clear sqrt ( egg ( y ) ) * z $
    for all y    let halfegg!_ ( y )     be sqrt ( y )     $
    for all y, z let halfegg!_ ( y ) * z be sqrt ( y ) * z $
    wrk1 takes wrk2 $
    for all y    clear halfegg!_ ( y )     $
    for all y, z clear halfegg!_ ( y ) * z $
    return wrk1 $
    end $
%
%
for all x such that part ( dummy1!_ * dummy2!_* x, 0 ) is minus 
  let egg ( x ) be egg ( x / ( - 1 ) ) * ( - 1 ) ** p!_ $
%
for all x, y let egg ( x * y ) be egg ( x ) * egg ( y ) $
%
for all x such that fixp ( den x ) and den x > 1  
  let egg ( x ) be egg ( num x ) / egg ( den x ) $ 
%
for all n such that numberp n and n > 0 let egg ( n ) be n $
%
% A function to facilitate substitutions of angular momentum variables.
% And also can be used simply to make the expressions neat.
%
% Note ( 1 ) This function behaves the same as the reduce function
%            "sub" except that the phase of the expression is
%            properly evaluated.
%
%      ( 2 ) This function does not cause any change on the input
%            expression to be evaluated, the second argument of
%            eval: "y".
%
%
operator eval $
%
for all x, y let eval ( x, y ) be
  begin scalar wrk1, wrk2 $
    reduction!_control!_ ( 0 ) $
    clear p!_, q!_ $
    reduction!_control!_ ( 2 ) $
    work!_for!_meta!_mode!_of!_krdlta!_ takes y $
    meta!_mode!_of!_krdlta!_ takes active $
    wrk2 takes y $
    meta!_mode!_of!_krdlta!_ takes inactive $
    wrk1 takes  work!_for!_meta!_mode!_of!_krdlta!_ $
    reduction!_control!_ ( 3 ) $
    reduction!_control!_ ( 4 ) $
    wrk2 takes sub ( x, wrk1 ) $
    wrk1 takes sub ( p!_ = q!_ , wrk2 ) $
    if phase!_ is concrete then do
      begin $
        let sqrt ( - 1 ) ** q!_ be i $
        q!_ takes 1 $
        end $
    wrk2 takes wrk1 $
    if phase!_ is concrete then do
      begin $
        clear q!_ $
        clear sqrt ( - 1 ) ** q!_ $
        end $
    wrk1 takes sub ( q!_ = p!_, wrk2 ) $
    reduction!_control!_ ( 8 ) $
    return wrk1 $
    end $
%
for all y let eval ( y ) be eval ( {}, y ) $
%
%
%
% End of Chapter 1. Collection of auxiliary operators and procedures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% Chapter 2. Clebsch-Gordan coefficients and related quantities        %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% The full path of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter20.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter20.r,v 1.1 92/04/16 20:37:43 koikef Locked $'/
%
%
%
% Wigner's formula for the Clebsch-Gordan coefficients
%
%      References:
%
%      1. E. P. Wigner, Gruppentheorie, Friedrich Vieweg und Sohn
%         Braunschweig, 1931.
%      2. E. P. Wigner, Group Theory, Academic Press, 1959.
%      3. E. P. Wigner, ( Morita Masato and Morita Reiko, Japanese
%         translation ), Gunron to ryoushirikigaku, yoshiokashoten.
%
%
%
%  Auxiliary functions:
%
%  Note: The phase variable p!_ must be cleared out before any
%        call of the following auxiliary functions.
%        Otherwise, the correctness of the output is not guaranteed.
%
%  Note: The following auxiliary functions are normally unnecessary
%        to be called directly by the user. No output for an input
%        afterwards the direct call of these auxiliary functions
%        will be guaranteed of their correctness.
%
%
%  Find a range of the summation in the Wigner's formula.
%
procedure findsumrangecg ( j1, j2, j3, m1, m2, m3 ) $
  begin $
    clear nus!_, nux!_, nun!_ $
    if fixp ( j3 - j1 + j2 ) then nus!_ takes j3 - j1 + j2 $
    if fixp ( j3 + j1 - j2 ) then nus!_ takes j3 + j1 - j2 $
    if fixp ( j3 + m3      ) then nus!_ takes j3 + m3      $
    if fixp ( j3 - m3      ) then nus!_ takes j3 - m3      $
    if fixp nus!_ then do
      begin $
        if fixp ( j3 - j1 + j2 ) and j3 - j1 + j2 <= nus!_ then do
          begin $
            nux!_ takes j3 - j1 + j2  $
            nun!_ takes 0             $
            nus!_ takes nux!_ - nun!_ $
            end $
        if fixp ( j3 + j1 - j2 ) and j3 + j1 - j2 <= nus!_ then do
          begin $
            nux!_ takes j3 + m3       $
            nun!_ takes j2 - j1 + m3  $
            nus!_ takes nux!_ - nun!_ $
            end $
        if fixp ( j3 + m3      ) and j3 + m3      <= nus!_ then do
          begin $
            nux!_ takes j3 + m3       $
            nun!_ takes 0             $
            nus!_ takes nux!_ - nun!_ $
            end $
        if fixp ( j3 - m3      ) and j3 - m3      <= nus!_ then do
          begin $
            nux!_ takes j3 - j1 + j2  $
            nun!_ takes j2 - j1 + m3  $
            nus!_ takes nux!_ - nun!_ $
            end $
        findsumrange takes successful $
        end
                  else do
      begin $
        clear nun!_, nux!_, nus!_ $
        findsumrange takes unsuccessful $
        end $
    return findsumrange $
    end $
%
%
%  The nu-sum of
%
%
%     nu+j2+m2    (j3+j1-j2)!(j3-j1+j2)!
%  (-)         * ------------------------
%                  (nu)!(j3-j1+j2-nu)!
%
%
%                    (j3+m3)!(j3-m3)!
%              * -------------------------
%                (j3+m3-nu)!(nu+j1-j2-m3)!
%
%
%                    (j2+j3+m1-nu)!(j1-m1+nu)!
%              * ---------------------------------
%                 (j2+j3+m1-nux!_)!(j1-m1+nun!_)!
%
%
%  over all the possible nu's that give non-negative
%  arguments for the factorials:
%
procedure cgaux!#!_02 ( j1, j2, j3, m1, m2, m3 ) $
  begin scalar nu, wrk, wrka, wrkb, wrk1, wrk2, wrk3, wrk4, cgauxwk $
    nu takes  nun!_ $
    p!_ takes 1 $
    wrk  takes 0 $
    while     nux!_ - nu >= 0
    do
      begin $
        on exp $
        wrk1 takes   ftcdiv!_ ( j3 + j1 - j2 , nu )
                   * ftcdiv!_ ( j3 - j1 + j2 , j3 - j1 + j2 - nu )
                   $
        wrk2 takes   ftcdiv!_ ( j3 + m3 , j3 + m3 - nu )
                   * ftcdiv!_ ( j3 - m3 , nu + j1 -j2 -m3 )
                   $
        wrk3 takes   ftcdiv!_ 
                       ( j2 + j3 + m1 - nu , j2 + j3 + m1 - nux!_ )
                   * ftcdiv!_ 
                       ( j1 - m1 + nu , j1 - m1 + nun!_ )
                   $
        wrk4 takes   wrk1 * wrk2 * wrk3
                   $
        wrk  takes   wrk  + ( - 1 ) ** ( ( nu - nun!_ ) * p!_ ) * wrk4
                   $
        nu takes nu + 1 $
      end $
    clear p!_ $
    if wrk isnot 0 and part ( dummy1!_ * dummy2!_ * wrk , 0 ) is minus
      then wrk takes wrk / ( - 1 ) * ( - 1 ) ** ( p!_ ) $
    cgauxwk takes ( - 1 ) ** ( p!_ * ( j2 + m2 + nun!_ ) ) * wrk $
    return cgauxwk $
    end $
%
%
%                  (2j3+1)                        (j1+j2-j3)!
%    ---------------------------------------- * ---------------
%     (j3+m3)!(j3-m3)!(j3+j1-j2)!(j3-j1+j2)!     (j1+j2+j3+1)!
%
%
%     (j2+j3+m1-nux!_)!(j2+j3+m1-nux!_)!(j1-m1+nun!_)!(j1-m1+nun!_)!
%  * ----------------------------------------------------------------
%                    (j1+m1)!(j2-m2)!(j1-m1)!(j2+m2)!
%
%
procedure cgaux!#!_01 ( j1, j2, j3, m1, m2, m3 ) $
  begin scalar wrk0,
               wrk1, wrk11, wrk12, wrk13,wrk14,wrk15,wrk16,
               wrk2,
               wrk3, wrk31, wrk32, wrk33,wrk34,wrk35,wrk36,
               wrk4 $
        on exp $
        wrk11 takes tactrl!_ ( j3 + m3 )      $
        wrk12 takes tactrl!_ ( j3 - m3 )      $
        wrk13 takes tactrl!_ ( j3 + j1 - j2 ) $
        wrk14 takes tactrl!_ ( j3 - j1 + j2 ) $
      wrk15 takes wrk11 * wrk12 $
      wrk16 takes wrk13 * wrk14 $
    wrk1 takes egg ( 2 * j3 + 1 ) * wrk15 * wrk16 $
    wrk2 takes  ftcdiv!_ ( j1 + j2 - j3 , j1 + j2 + j3 + 1 ) $
        wrk31 takes  ftcdiv!_ ( j2 + j3 + m1 - nux!_ , j1 + m1 ) $
        wrk32 takes  ftcdiv!_ ( j2 + j3 + m1 - nux!_ , j2 - m2 ) $
        wrk33 takes  ftcdiv!_ ( j1 - m1 + nun!_ , j1 - m1 )      $
        wrk34 takes  ftcdiv!_ ( j1 - m1 + nun!_ , j2 + m2 )      $
      wrk35 takes wrk31 * wrk32 $
      wrk36 takes wrk33 * wrk34 $
    wrk3 takes wrk35 * wrk36 $
    wrk4 takes wrk1 * wrk2 $
    wrk0 takes wrk3 * wrk4 $
    return wrk0 $
    end $
%
%
%  Square of Wigner's formula
%
procedure cgaux!#!_00 ( j1, j2, j3, m1, m2, m3 ) $
  begin scalar wrk0, wrk1, wrk2, wrk3 $
    if    not numberp j1 or not numberp j2 or not numberp j3
       or not numberp m1 or not numberp m2 or not numberp m3 then do
      begin $
        reduction!_control!_ ( 1 ) $
        on exp $
        wrk1 takes cgaux!#!_01 ( j1, j2, j3, m1, m2, m3 ) $
        on exp $
        wrk2 takes cgaux!#!_02 ( j1, j2, j3, m1, m2, m3 ) $
        on exp $
        wrk3 takes wrk1 * wrk2 $
        wrk0 takes wrk2 * wrk3 $
        reduction!_control!_ ( 2 ) $
        wrk1 takes wrk0 $
        reduction!_control!_ ( 3 ) $
        on exp $
        wrk2 takes wrk1 $
        end
                                                             else do
      begin $
        p!_ takes 1 $
        wrk1 takes cgaux!#!_01 ( j1, j2, j3, m1, m2, m3 ) $
        wrk2 takes cgaux!#!_02 ( j1, j2, j3, m1, m2, m3 ) $
        clear p!_ $
        if     wrk2 isnot 0 
           and part ( wrk2 * dummy1!_ * dummy2!_, 0 ) is minus then
           wrk2 takes wrk2 / ( - 1 ) * ( - 1 ) ** p!_ $
        wrk3 takes wrk2 ** 2 $
        wrk0 takes wrk1 * wrk3 $
        wrk2 takes wrk0 $
        end $
    return wrk2 $
    end $
%
%
%  Check (1) if the triangular conditions are fulfilled, and
%        (2) if ji and mi, where i = 1, 2, 3, are a pair of
%            legal angular momentum quantum numbers.
%
procedure trangcondcg ( j1, j2, j3, m1, m2, m3 ) $
  begin scalar flag $
    if (
            fixp ( j1 - m1 ) and ( j1 - m1 ) < 0
         or fixp ( j1 + m1 ) and ( j1 + m1 ) < 0
      or
            fixp ( j2 - m2 ) and ( j2 - m2 ) < 0
         or fixp ( j2 + m2 ) and ( j2 + m2 ) < 0
      or
            fixp ( j3 - m3 ) and ( j3 - m3 ) < 0
         or fixp ( j3 + m3 ) and ( j3 + m3 ) < 0
      or
            fixp ( j1 + j2 - j3 ) and ( j1 + j2 - j3 ) < 0
      or
            fixp ( j1 - j2 + j3 ) and ( j1 - j2 + j3 ) < 0
      or
            fixp (-j1 + j2 + j3 ) and (-j1 + j2 + j3 ) < 0
      or
            fixp ( ( 2 * ( j1 + j2 + j3 ) + 1 ) / 2 )
      or
            fixp ( ( 2 * ( j1 - m1 ) + 1 ) / 2 )
      or
            fixp ( ( 2 * ( j2 - m2 ) + 1 ) / 2 )
      or
            fixp ( ( 2 * ( j3 - m3 ) + 1 ) / 2 )
       )
        then flag takes denied
        else flag takes affirmed
      $
    triangular!_cond!_check!_mode!_ takes Clebsch $
    Clebsch!_trang!_cond!_1!_save!_ takes   j1 + j2 - j3 $
    Clebsch!_trang!_cond!_2!_save!_ takes   j1 - j2 + j3 $
    Clebsch!_trang!_cond!_3!_save!_ takes - j1 + j2 + j3 $
    Clebsch!_trang!_cond!_4!_save!_ takes j1 - m1 $
    Clebsch!_trang!_cond!_5!_save!_ takes j2 - m2 $
    Clebsch!_trang!_cond!_6!_save!_ takes j3 - m3 $
    Clebsch!_trang!_cond!_7!_save!_ takes j1 + m1 $
    Clebsch!_trang!_cond!_8!_save!_ takes j2 + m2 $
    Clebsch!_trang!_cond!_9!_save!_ takes j3 + m3 $
    Clebsch!_trang!_cond!_a!_save!_ takes j1  $
    Clebsch!_trang!_cond!_b!_save!_ takes j2  $
    Clebsch!_trang!_cond!_c!_save!_ takes j3  $
    return flag $
    end $
%
%
%  Find one of the best conditions for calculating the
%  C. G. coefficient with a given set of parameters, and
%  perform the calculation.
%
%
procedure cgc!_aux!_proc!_ 
  ( cgc!_aux!_index!_, j1, j2, j3, m1, m2, m3 ) $
  begin scalar wrk0, wrk1, wrk2, wrk3, flag $
    cgcoefficient!_indexw!_ takes undetermined $
    cj1!_ takes j1 $ cj2!_ takes j2 $ cj3!_ takes j3 $
    cm1!_ takes m1 $ cm2!_ takes m2 $ cm3!_ takes m3 $
    wrk0 takes krdlta ( m3, m1 + m2 ) $
    flag takes trangcondcg ( j1, j2, j3, m1, m2, m3 ) $
    if flag is affirmed and wrk0 isnot 0 then do
      begin $
        ww!_x!_ takes wrk0 $
        cj1!_ takes make!_window!_for!_x!_ ( j1 ) $
        cj2!_ takes make!_window!_for!_x!_ ( j2 ) $
        cj3!_ takes make!_window!_for!_x!_ ( j3 ) $
        cm1!_ takes make!_window!_for!_x!_ ( m1 ) $
        cm2!_ takes make!_window!_for!_x!_ ( m2 ) $
        cm3!_ takes cm1!_ + cm2!_ $
        wrk0 takes ww!_x!_ $
        clear nus!_m!_, nux!_m!_, nun!_m!_ $
        cgcoefficient!_index!_ takes cgc!_aux!_proc!_aux!_
          ( 1, cj1!_, cj2!_, cj3!_, cm1!_, cm2!_, cm3!_ ) $
        cgcoefficient!_index!_ takes cgc!_aux!_proc!_aux!_
          ( 2, cj2!_, cj3!_, cj1!_,-cm2!_, cm3!_, cm1!_ ) $
        cgcoefficient!_index!_ takes cgc!_aux!_proc!_aux!_
          ( 3, cj3!_, cj1!_, cj2!_, cm3!_,-cm1!_, cm2!_ ) $
        nus!_ takes nus!_m!_ $ 
        nux!_ takes nux!_m!_ $ 
        nun!_ takes nun!_m!_ $
        end
                                                          else do
      begin $
        cgcoefficient!_index!_ takes 0 $
        end $
    if cgcoefficient!_index!_ isnot undetermined then do
      begin $
        clear p!_, q!_ $
        reduction!_control!_ ( 0 ) $
        wrk2 takes c2!_ ( cj1!_, cj2!_, cj3!_, cm1!_, cm2!_, cm3!_ ) $
        on exp $
        if cgc!_aux!_index!_ is 1 then do
          begin $
            wrk1 takes sqrt ( num wrk2 ) $
            wrk3 takes sqrt ( den wrk2 ) $
            wrk2 takes hatch!_half!_eggs!_ ( wrk3 ) $
            wrk3 takes hatch!_half!_eggs!_ ( wrk1 ) $
            wrk1 takes wrk3 / wrk2 $
            end $
        if cgc!_aux!_index!_ is 2 then do
          begin $
            wrk1 takes wrk2 $
            end $
        if cgc!_aux!_index!_ is 3 then do
          begin $
            wrk3 takes      ( wrk2
              / ( - 1 ) ** ( p!_ * ( cj1!_ - cj2!_ + cm3!_ ) * 2 )
              / egg ( 2 * cj3!_ + 1 ) ) $
            wrk1 takes sqrt ( num wrk3 ) $
            wrk2 takes sqrt ( den wrk3 ) $
            wrk3 takes hatch!_half!_eggs!_ ( wrk2 ) $
            wrk2 takes hatch!_half!_eggs!_ ( wrk1 ) $
            wrk1 takes wrk2 / wrk3 $
            end $
        if cgc!_aux!_index!_ is 4 then do
          begin $
            wrk1 takes        wrk2
              / ( - 1 ) ** ( p!_ * ( cj1!_ - cj2!_ + cm3!_ ) * 2 )
              / egg ( 2 * cj3!_ + 1 )   $
            end $
        if cgc!_aux!_index!_ is 5 then do
          begin $
            wrk3 takes      ( wrk2
              / ( - 1 ) ** ( p!_ * ( ( cj3!_  + cm3!_ ) * 2 ) )
              / egg ( 2 * cj3!_ + 1 ) ) $
            wrk1 takes sqrt ( num wrk3 ) $
            wrk2 takes sqrt ( den wrk3 ) $
            wrk3 takes hatch!_half!_eggs!_ ( wrk2 ) $
            wrk2 takes hatch!_half!_eggs!_ ( wrk1 ) $
            wrk1 takes wrk2 / wrk3 $
            end $
        if cgc!_aux!_index!_ is 6 then do
          begin $
            wrk1 takes        wrk2
              / ( - 1 ) ** ( p!_ * ( ( cj3!_  + cm3!_ ) * 2 ) )
              / egg ( 2 * cj3!_ + 1 )   $
            end $
        if cgc!_aux!_index!_ is 7 then do
          begin $
            wrk3 takes      ( wrk2
              * ( - 1 ) ** ( p!_ * ( ( cj1!_ + cj2!_ + cj3!_ ) * 2 ) )
              / ( - 1 ) ** ( p!_ * ( ( cj1!_ - cj2!_ + cm3!_ ) * 2 ) )
              / egg ( 2 * cj3!_ + 1 ) ) $
            wrk1 takes sqrt ( num wrk3 ) $
            wrk2 takes sqrt ( den wrk3 ) $
            wrk3 takes hatch!_half!_eggs!_ ( wrk2 ) $
            wrk2 takes hatch!_half!_eggs!_ ( wrk1 ) $
            wrk1 takes wrk2 / wrk3 $
            end $
        if cgc!_aux!_index!_ is 8 then do
          begin $
            wrk1 takes        wrk2
              * ( - 1 ) ** ( p!_ * ( ( cj1!_ + cj2!_ + cj3!_ ) * 2 ) )
              / ( - 1 ) ** ( p!_ * ( ( cj1!_ - cj2!_ + cm3!_ ) * 2 ) )
              / egg ( 2 * cj3!_ + 1 )   $
            end $
        if phase!_ is concrete then p!_ takes 1 $
        wrk2 takes wrk1 * wrk0 $
        clear p!_, q!_ $
        reduction!_control!_ ( 8 ) $
        end
                                              else do
      begin $
        wrk2 takes program!_error $
        end $
      off exp $
      wrk3 takes hatcheggs ( wrk2 ) $
      off exp $
    return wrk3 $
    end $
%
%
procedure cgc!_aux!_proc!_aux!_ 
  ( index!_no!_, j1, j2, j3, m1, m2, m3 ) $
  begin scalar flag $
    flag takes findsumrangecg ( j1, j2, j3, m1, m2, m3 ) $
    if flag is successful then do
      begin $
        if    not fixp nus!_m!_ and fixp nus!_
           or fixp nus!_m!_ and nus!_m!_ > nus!_ then do
          begin $
            nus!_m!_ takes nus!_ $
            nun!_m!_ takes nun!_ $
            nux!_m!_ takes nux!_ $
            cgcoefficient!_indexw!_ takes index!_no!_ $
            end $
      end $
    return cgcoefficient!_indexw!_ $
    end $
%
%
%  An auxiliary operator c2!_ to represent the Clebsch-Gordan
%  coefficients.
%
operator c2!_ $
%
for all j1, j2, j3, m1, m2, m3 such that
  cgcoefficient!_index!_ is 2
  let c2!_ ( j1, j2, j3, m1, m2, m3 ) be
    begin scalar wrk1, wrk2 $
      wrk1 takes   ( - 1 ) ** ( p!_ * ( ( j2 + m2 ) * 2 ) )
              * ( 2 * j3 + 1 ) / ( 2 * j1 + 1 )
              * cgaux!#!_00 ( j2, j3, j1, - m2, m3, m1 ) $
      wrk2 takes wrk1 $
      return wrk2 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that
  cgcoefficient!_index!_ is 3
  let c2!_ ( j1, j2, j3, m1, m2, m3 ) be
    begin scalar wrk1, wrk2 $
      wrk1 takes   ( - 1 ) ** ( p!_ * ( ( j1 - m1 )  * 2 ) )
              * ( 2 * j3 + 1 ) / ( 2 * j2 + 1 )
              * cgaux!#!_00 ( j3, j1, j2, m3, -m1, m2 ) $
      wrk2 takes wrk1 $
      return wrk2 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that
  cgcoefficient!_index!_ is 1
  let c2!_ ( j1, j2, j3, m1, m2, m3 ) be
    begin scalar wrk1, wrk2 $
      wrk1 takes cgaux!#!_00 ( j1, j2, j3, m1, m2, m3 ) $
      wrk2 takes wrk1 $
      return wrk2 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that
  cgcoefficient!_index!_ is 0
  let c2!_ ( j1, j2, j3, m1, m2, m3 ) be 0 $
%
%
%  End of  Auxiliary functions.
%
%
%
%  Note: In the following entries, the phase variable p!_
%        is left abstract on exit if phase!_ is set abstract,
%        and the p!_ is set 1 on exit if phase!_ is set concrete.
%        This feature may be controlled by using a function 
%        "phasemode".
%
%        Input
%
%          phasemode concrete ;
%
%        to get an expression withOUT an undetermined phase variable 
%        p!_. And input
%
%          phasemode abstract ;
%
%        to get an expression with    an undetermined phase variable 
%        p!_. The default is "phase concrete ;".
%
%
%  2.1 Clebsch-Gordan coefficients in terms of the Rose notation.
%
%    References:
%
%    1. M. E. Rose, Elementary Theory of Angular Momentum,
%       John Wiley and Sons Inc., N. Y., 1957
%    2. M. E. Rose, ( Yamanouchi Takahiko and Morita Masato,
%       Japanese translation), Kakuundouryouno kisoriron,
%       Misuzu shobou, 1971.
%
%
%  Rose's c - coefficients :
%
operator c $
%
for all j1, j2, j3, m1, m2 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 )
     or fixp ( j3 + m1 + m2 ) or fixp ( j3 - m1 - m2 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let c ( j1, j2, j3, m1, m2 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 1, j1, j2, j3, m1, m2, m1 + m2 ) $
      return wrk0 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 )      or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let c ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 1, j1, j2, j3, m1, m2, m3 ) $
      return wrk0 $
      end $
%
%  Square of the Rose's c - coefficients :
%
operator c2 $
%
for all j1, j2, j3, m1, m2 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 )
     or fixp ( j3 + m1 + m2 ) or fixp ( j3 - m1 - m2 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let c2 ( j1, j2, j3, m1, m2 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 2, j1, j2, j3, m1, m2, m1 + m2 ) $
      return wrk0 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 )      or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let c2 ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 2, j1, j2, j3, m1, m2, m3 ) $
      return wrk0 $
      end $
%
%
% 2.2 The Wigner 3-j symbols
%
%    References:
%
%    1. E. P. Wigner, in L. C. Biedenharn and H. van Dam, eds.,
%       Quantum theory of Angular momentum ( Academic Press
%       New York, 1965 ).
%
%    2. A. R. Edmonds, Angular momentum in quantum mechanics
%       ( Princeton University Press, Princeton, N. J.,1960),
%       2'nd ed.
%    3. R. D. Cowan, The theory of atomic structure and spectra
%       ( University of California Press, Berkeley, Los Angels,
%       London , 1981 ) in Chap. 5.
%
%
% Wigner's 3j-symbols:
%
operator w3j $
%
%                                      j1, j2, j3
%  w3j ( j1, j2, j3, m1, m2, m3 ) =  (            )
%                                      m1, m2, m3
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 )      or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let w3j ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 3, j1, j2, j3, m1, m2, - m3 ) $
      return wrk0 $
      end $
%
%
% Square of Wigner's 3j-symbols:
%
operator w3j2 $
%
%                                         j1, j2, j3    2
%  w3j2 ( j1, j2, j3, m1, m2, m3 ) =  [ (            ) ]
%                                         m1, m2, m3
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 ) or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 ) or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 ) or fixp ( j2 - m2 ) )
  let w3j2 ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 4, j1, j2, j3, m1, m2, - m3 ) $
      return wrk0 $
      end $
%
%
% 2.3 Racah's V-coefficients and the V-bar coefficients
%     by Fano and Racah
%
%    References:
%
%    1. G. Racah, II, Phys. Rev. 62, 438 (1942).
%
%    2. U. Fano and G. Racah, Irreducible tensorial sets
%       ( Academic Press, New York, 1959 ).
%
%
% Racah's V-coefficients:
%
operator v $
%
%  v ( j1, j2, j3, m1, m2, m3 )
%
%       -j3+m3      -1/2
%   = (-)     (2j3+1)   c ( j1, j2, j3, m1, m2, -m3 )
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 )      or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let v ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 5, j1, j2, j3, m1, m2, - m3 ) $
      return wrk0 $
      end $
%
% Square of Racah's V-coefficients:
%
operator v2 $
%
%        2
%  v2 = v
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 )      or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let v2 ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 6, j1, j2, j3, m1, m2, - m3 ) $
      return wrk0 $
      end $
%
% Racah's V-bar coefficients
%
operator vbar $
%
%            2(j2+j3)       j1+j2+j3
%  vbar = (-)        v = (-)        w3j
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 )      or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let vbar ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 7, j1, j2, j3, m1, m2, - m3 ) $
      return wrk0 $
      end $
%
% Square of Racah's V-bar coefficients:
%
operator vbar2 $
%
%              2
%  vbar2 = vbar
%
for all j1, j2, j3, m1, m2, m3 such that
        symbolic!_mode!_ is inactive
  and   cgc!_red!_mode!_00 is active
  and ( fixp ( j3 - j1 + j2 ) or fixp ( j3 + j1 - j2 ) 
     or fixp ( j3 + m3 )      or fixp ( j3 - m3 )
     or fixp ( j1 - j2 + j3 ) or fixp ( j1 + j2 - j3 ) 
     or fixp ( j1 + m1 )      or fixp ( j1 - m1 )
     or fixp ( j2 - j3 + j1 ) or fixp ( j2 + j3 - j1 ) 
     or fixp ( j2 + m2 )      or fixp ( j2 - m2 ) )
  let vbar2 ( j1, j2, j3, m1, m2, m3 ) be
    begin wrk0 $
      wrk0 takes cgc!_aux!_proc!_ ( 8, j1, j2, j3, m1, m2, - m3 ) $
      return wrk0 $
      end $
%
%
% 2.4 Transformation between C. G. coefficients and Wigner's 3j-symbols
%
%
% Notes: (0)   The following features will be effective only
%            when the variable symbolic!_mode!_ is being set active.
%            In this mode the program does not reduce
%            angular_momentum coefficients into concrete form.
%            This mode can be set by calling a procedure evaluation, as
%
%              evaluation deactivate ;
%
%            The default is "activate", with which the reductions
%            will always be performed for the coefficients and symbols,
%            if possible.
%
%        (1)   When the variable c!_to!_w3j!_ is set "active"
%            and the variable w3j!_to!_c!_ is set "inactive"
%            at the call of c or c2, they are transformed into
%            the corresponding w3j or w3j2 prior to the
%            reduction of the coefficients. This mode can be
%            set by calling a procedure notation, as
%
%              notation Wigner ;
%
%            The default is "Asitis", with which no transformation
%            will be performed between the C. G.'s and 3j-symbols.
%
%        (2)   When the variable c!_to!_w3j!_ is set "inactive"
%            and the variable w3j!_to!_c!_ is set "active"
%            at the call of  w3j or w3j2, they are transformed into
%            the corresponding c or c2 prior to the
%            reduction of the coefficients. This mode can be
%            set by calling a procedure notation, as
%
%              notation Rose ;
%
%            The default is "Asitis", with which NO transformation
%            will be performed between the C. G.'s and 3j-symbols.
%
%        (3)   When both the variables  c!_to!_w3j!_ and w3j!_to!_c!_
%            are set "inactive", NO transformation will be performed.
%            This mode is the default and can be reset by calling a
%            procedure notation, as
%
%              notation Asitis ;
%
%        (4)   When both the variables  c!_to!_w3j!_ and w3j!_to!_c!_
%            are set "active", NO transformation will be performed
%            at all.
%
%
for all j1, j2, j3, m1, m2, m3 such that c!_to!_w3j!_ is   active
                                     and w3j!_to!_c!_ is inactive
                                     and symbolic!_mode!_ is active
  let c ( j1, j2, j3, m1, m2, m3 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                     ( - 1 ) ** ( p!_ * ( j1 - j2 + m3 ) )
                   * sqrt ( 2 * j3 + 1 )
                   * w3j ( j1, j2, j3, m1, m2, -m3 ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that c!_to!_w3j!_ is inactive
                                     and w3j!_to!_c!_ is   active
                                     and symbolic!_mode!_ is active
  let w3j ( j1, j2, j3, m1, m2, m3 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                     ( - 1 ) ** ( p!_ * ( - j1 + j2 + m3 ) )
                   / sqrt ( 2 * j3 + 1 )
                   * c ( j1, j2, j3, m1, m2, -m3 ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that c!_to!_w3j!_ is   active
                                     and w3j!_to!_c!_ is inactive
                                     and symbolic!_mode!_ is active
  let c2 ( j1, j2, j3, m1, m2, m3 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                     ( - 1 ) ** ( p!_ * ( ( j1 - j2 + m3 ) * 2 ) )
                   * ( 2 * j3 + 1 )
                   * w3j2 ( j1, j2, j3, m1, m2, -m3 ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
for all j1, j2, j3, m1, m2, m3 such that c!_to!_w3j!_ is inactive
                                     and w3j!_to!_c!_ is   active
                                     and symbolic!_mode!_ is active
  let w3j2 ( j1, j2, j3, m1, m2, m3 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                     ( - 1 ) ** ( p!_ * ( ( - j1 + j2 + m3 ) * 2 ) )
                   / ( 2 * j3 + 1 )
                   * c2 ( j1, j2, j3, m1, m2, -m3 ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
%
% End of Chapter 2. Clebsch-Gordan coefficients and related quantities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% Chapter 3. Useful quantities related to the C. G. coefficients       %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% The full path name of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter30.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter30.r,v 1.1 92/04/16 20:37:46 koikef Locked $'/
%
%
% 3.1. Transformation coefficient of the hydrogenic atomic
%      state wavefunctions between the polar coordinate
%      representation and the parabolic coordinate representation.
%
%    References:
%
%    1. J. W. B. Huges, Proc. Phys. Soc., 91, 810 (1967), "Stark
%       states and O(4) symmetry of hydrogenic atoms".
%
%    2. David Park, Z. Phys. 159 155 (1960), "Relation between
%       the Parabolic and Special Eigenfunctions of Hydrogen".
%
%      trpopa ( n, q, l, m ) = < n q m | l m n >
%
%         (n+q-m-1)/2   n-1   n-1     m+q   m-q
%    = (-)           C(-----,-----,l,-----,-----,m)
%                        2     2       2     2
%
operator trpopa $
for all n, q, l, m such that fixp l and fixp m
  let trpopa ( n, q, l, m ) be
    begin scalar cwk, cwk1, cwk2, cwk3, stat, psave $
      psave takes phase!_ $
      clear p!_ $
      clear a, b $
      cwk takes   c2 ( a, a, l, m + b, - b, m ) $
      cwk1 takes cwk * ( - 1 ) ** ( p!_ * ( 2 * a + 2 * b ) ) $
      cwk2 takes sqrt ( cwk1 ) $
      a takes ( n - 1 ) / 2 $
      b takes ( q - m ) / 2 $
      cwk3 takes  cwk2 $
      if psave is concrete then p!_ takes 1 $
      stat takes cwk3 $
      clear p!_ $
      phase!_ takes psave $
      return stat $
      end $
for all n, q, l, m such that not ( fixp l and fixp m )
  let trpopa ( n, q, l, m ) be
    begin scalar a, b, cwk1, stat $
      clear p!_ $
      a takes ( n - 1 ) / 2 $
      b takes ( q - m ) / 2 $
      cwk1 takes      (
                       ( - 1 ) ** ( p!_ * ( ( 2 * a + 2 * b ) / 2 ) )
                     * c ( a, a, l, m + b, - b, m )
                    ) $
      if phase!_ is concrete then p!_ takes 1 $
      stat takes cwk1 $
      clear p!_ $
      return stat $
      end $
%
% Generate a fortran program that gives the trpopa.
%
procedure trpopafort ( lmax ) $
begin scalar l, m, stat, nmax $
nmax takes lmax + 1 $
phase!_ takes concrete $
on fort $
fortwidth!* takes 60 $
off period $
out "trpopa.fort" $
%write "c    To make a program with a floating " $
%write "c    point data table, " $
%write "c    run this program" $
%write "c    after replacing" $
%write "c" $
%write "c    (-1.)**N by (-1)**NN" $
%write "c" $
%write "c    .) by .D0)  ." $
%write "c    Please do them in the order stated." $
%write "c" $
 write "c  CUT HERE." $
 write "c>--------------------------------------<" $
%write "      program main" $
%write "      implicit real *8 (a-h,o-z)" $
%write "      integer n,l,m,q,qmin,qmax" $
%write "      character *1 a(80)" $
%write "      data nmax / ",nmax,"  /" $
%write "      open ( 8, file = 'trpopa.template_1' ) " $
%write "      open ( 9, file = 'trpopa.template_2' ) " $
%write "      do 5 i = 1,1000" $
%write "      read(8,101,end=7) ( a(j), j = 1, 80 )" $
%write "      write(6,101)      ( a(j), j = 1, 80 )" $
%write "    5 continue" $
%write "    7 continue" $
%write "      icount = 0" $
%write "      do 10 n = 1, nmax" $
%write "        do 10 m = 0, n-1" $
%write "          do 10 l = m, n-1" $
%write "            qmin = -( n - abs(m) - 1 )" $
%write "            qmax =  ( n - abs(m) - 1 )" $
%write "            do 10 q = qmin, qmax, 2" $
%write "              if ( q .ge. 0 ) then" $
%write "              icount = icount + 1" $
%write "              x = statel ( n, l, m, q )" $
%write "              write ( 6, 102 ) icount, x" $
%write "              end if " $
%write "   10 continue" $
%write "      do 15 i = 1,1000" $
%write "      read(9,101,end=17) ( a(j), j = 1, 80 )" $
%write "      write(6,101)       ( a(j), j = 1, 80 )" $
%write "   15 continue" $
%write "   17 continue" $
%write "      stop" $
%write "  101 format(80a1)" $
%write "  102 format(6x,'statmx( ',i6,' ) =',1pd24.16)" $
%write "      end" $
write "      function statel ( nn, l, m, qq )" $
write "c" $
write "c    This function gives the elements" $
write "c    of the transformation" $
write "c    matrix between the atomic " $
write "c    and stark states" $
write "c    for all possible combination" $
write "c    of the states with" $
write "c    the principal quantum number " $
write "c    n =< ", nmax,"." $
write "c" $
write "      integer nn, l, m, qq" $
write "      real *8 n, q" $
write "      real *8 statel" $
write "c" $
write "      n = nn" $
write "      q = qq" $
  clear n, q $
  for m takes 0 : lmax do
  begin $
    write "       if ( m .eq. ", m, ") then" $
    for l takes m : lmax do
    begin $
      write "        if ( l .eq. ", l, ") then" $
      stat takes trpopa ( n, q, l, m ) $
      on period $
      write "          statel = ", stat $
      off period $
      write "          return" $
      write "         end if" $
      clear p!_ $
    end $
    write "       end if" $
  end $
write "      end" $
write "c>--------------------------------------<" $
write "c  CUT HERE." $
write "c" $
shut "trpopa.fort" $
end $
%
%
% End of Chapter 3. Useful quantities related to the C. G. coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% Chapter 4. Racah coefficients and Wigner's 6j-Symbols                %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% The full path name of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter40.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter40.r,v 1.1 92/04/16 20:37:50 koikef Locked $'/
%
%
% 4.1. Auxiliary operators and procedures
%
%
% Note: A phase variable p!_ must be cleared out before any call of
%       the following auxiliary procedures. Otherwise the correctness
%       of the output may not be guaranteed.
%
% Note: A direct call of the following auxiliary functions from the
%       top level is normally unnecessary. No outputs for an input
%       afterwards this direct call will be guaranteed of
%       their correctness.
%
%
% Modified triangular coefficients:
%
% deltar (a,b,c) = [{(a-b+c)!(-a+b+c)!}/{(a+b-c)!(a+b+c+1)!}] ** 1/2
%
% deltar2(a,b,c) = [{(a-b+c)!(-a+b+c)!}/{(a+b-c)!(a+b+c+1)!}]
%
procedure deltar ( j1, j2, j3 ) $
  begin scalar wrk1, wrk2, wrk3 $
    wrk1 takes deltar2 ( j1, j2, j3 ) $
    wrk2 takes sqrt ( wrk1 ) $
    wrk3 takes wrk2 $
    return wrk3 $
    end $
%
procedure deltar2 ( j1, j2, j3 ) $
  begin scalar wrk1, wrk2, wrk3, wrk4, wrk5, wrk6 $
    if (     ( numberp ( j1 + j2 ) and numberp j3
               and   j1 + j2 < j3 )
         or  ( numberp ( j2 + j3 ) and numberp j1
               and   j2 + j3 < j1 )
         or  ( numberp ( j3 + j1 ) and numberp j2
               and   j3 + j1 < j2 ) )
         then do
      begin $
        wrk2 takes 0 $
        end
         else do
      begin $
        on exp $
        wrk3 takes factrl!_ (   j1 - j2 + j3     ) $
        wrk4 takes factrl!_ ( - j1 + j2 + j3     ) $
        wrk5 takes tactrl!_ (   j1 + j2 - j3     ) $
        wrk6 takes tactrl!_ (   j1 + j2 + j3 + 1 ) $
        wrk1 takes wrk3 * wrk5 $
        wrk2 takes wrk4 * wrk6 $
        wrk3 takes wrk1 * wrk2 $
        wrk2 takes wrk3 $
        end $
    return wrk2 $
    end $
%
%
% Square of the product of modified triangular coefficients:
%
% waux!#00!_(abcdef) = deltar2(abe)deltar2(cde)deltar2(acf)deltar2(bdf)
%
procedure waux!#00!_ ( j1, j2, j3, j4, j5, j6 ) $
  begin scalar wrk0, wrk1, wrk2, wrk3, wrk4, wprod $
    on exp $
    wrk1 takes deltar2 ( j1, j2, j5 ) $
    if wrk1 is 0 then wprod takes 0
                 else do
      begin $
        wrk2 takes deltar2 ( j3, j4, j5 ) $
        if wrk2 is 0 then wprod takes 0
                     else do
          begin $
            wrk3 takes deltar2 ( j1, j3, j6 ) $
            if wrk3 is 0 then wprod takes 0
                         else do
              begin $
                wrk4 takes deltar2 ( j2, j4, j6 ) $
                if wrk4 is 0 then wprod takes 0
                             else do
                  begin $
                    wprod takes   wrk1 * wrk2 * wrk3 * wrk4 $
                    end $
                end $
            end $
        end $
    wrk0 takes wprod $
    return wrk0 $
    end $
%
%
% k-Sum of [(-)**k/k!]
%          *[(a+b+c+d+1-k)!]/[(e+f-b-c+k)!(e+f-a-d+k)!]
%          *[(a+b-e)!/(a+b-e-k)!]
%          *[(c+d-e)!/(c+d-e-k)!]
%          *[(a+c-f)!/(a+c-f-k)!]
%          *[(b+d-f)!/(b+d-f-k)!]
%
procedure waux!#01!_ ( j1, j2, j3, j4, j5, j6 ) $
  begin scalar wrk0, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, wsum $
    k!_ takes kpmin!_ $
    wsum takes 0 $
    while kpmax!_ - k!_ >= 0 do
      begin $
        on exp $
        wrk2 takes  ( - 1 ) ** ( p!_ * ( k!_ - kpmin!_ ) )
                    * tactrl!_ 
                        (                                      k!_ )
                    * factrl!_ 
                        (   j1 + j2 + j3 + j4            + 1 - k!_ )
                    * tactrl!_ 
                        (      - j2 - j3      + j5  + j6     + k!_ )
                    * tactrl!_ 
                        ( - j1           - j4 + j5  + j6     + k!_ )
                    $
        wrk3 takes    ftcdiv!_ ( j1 + j2 - j5 , j1 + j2 - j5 - k!_ ) $
        wrk4 takes    ftcdiv!_ ( j3 + j4 - j5 , j3 + j4 - j5 - k!_ ) $
        wrk5 takes    ftcdiv!_ ( j1 + j3 - j6 , j1 + j3 - j6 - k!_ ) $
        wrk6 takes    ftcdiv!_ ( j2 + j4 - j6 , j2 + j4 - j6 - k!_ ) $
        on exp $
        wsum takes    wsum + wrk2 * wrk3 * wrk4 * wrk5 * wrk6 $
        k!_ takes k!_ + 1 $
        end $
    wrk0 takes wsum * ( - 1 ) ** ( p!_ * kpmin!_ ) $
    return wrk0 $
    end $
%
%
% Find the most convenient range for k!_ summation in Racah's formula.
%
procedure findsumrangeracah ( j1, j2, j3, j4, j5, j6 ) $
  begin scalar flag $
    clear kpdel!_ $
    flag takes findsumrraux
   ( j1 + j2 + j3 + j4           + 1 ,                           0 ) $
    flag takes findsumrraux
   ( j1 + j2 + j3 + j4           + 1 ,      j2 + j3      - j5 - j6 ) $
    flag takes findsumrraux
   ( j1 + j2 + j3 + j4           + 1 , j1           + j4 - j5 - j6 ) $
    flag takes findsumrraux
   (           j3 + j4 - j5          ,                           0 ) $
    flag takes findsumrraux
   (           j3 + j4 - j5          ,      j2 + j3      - j5 - j6 ) $
    flag takes findsumrraux
   (           j3 + j4 - j5          , j1           + j4 - j5 - j6 ) $
    flag takes findsumrraux
   ( j1 + j2           - j5          ,                           0 ) $
    flag takes findsumrraux
   ( j1 + j2           - j5          ,      j2 + j3      - j5 - j6 ) $
    flag takes findsumrraux
   ( j1 + j2           - j5          , j1           + j4 - j5 - j6 ) $
    flag takes findsumrraux
   (      j2      + j4      - j6     ,                           0 ) $
    flag takes findsumrraux
   (      j2      + j4      - j6     ,      j2 + j3      - j5 - j6 ) $
    flag takes findsumrraux
   (      j2      + j4      - j6     , j1           + j4 - j5 - j6 ) $
    flag takes findsumrraux
   ( j1      + j3           - j6     ,                           0 ) $
    flag takes findsumrraux
   ( j1      + j3           - j6     ,      j2 + j3      - j5 - j6 ) $
    flag takes findsumrraux
   ( j1      + j3           - j6     , j1           + j4 - j5 - j6 ) $
    if fixp kpdel!_ then findsumrangeracah takes   successful
                    else findsumrangeracah takes unsuccessful $
    return findsumrangeracah $
    end $
%
procedure findsumrraux ( kp1!_, kp2!_ ) $
  begin flag $
    flag takes unsuccessful $
    if fixp ( kp1!_ - kp2!_ ) then do
      begin $
        if  not fixp kpdel!_ 
           or ( fixp kpdel!_ and kp1!_ - kp2!_ < kpdel!_ ) then do
          begin $
            kpdel!_ takes kp1!_ - kp2!_ $
            kpmax!_ takes kp1!_ $
            kpmin!_ takes kp2!_ $
            flag takes successful $
            end $
      end $
  return flag $
  end $
%
%
% Check if the triangular conditions are fulfilled.
% Save triangular information for later references.
%
% procedure tragcondrc returns
%
%           'affirmed' if the condition is fulfilled or
%                      if the condition can not be checked,
%
%           or
%
%           'denied'   if the condition is NOT fulfilled.
%
procedure tragcondrc ( j1, j2, j3, j4, j5, j6 ) $
  begin scalar flag $
    if
     (
            (        fixp    ( j1 + j2 - j5 ) and j1 + j2 - j5 >= 0
              or not numberp ( j1 + j2 - j5 ) )
       and
            (        fixp    ( j2 + j5 - j1 ) and j2 + j5 - j1 >= 0
              or not numberp ( j2 + j5 - j1 ) )
       and
            (        fixp    ( j5 + j1 - j2 ) and j5 + j1 - j2 >= 0
              or not numberp ( j5 + j1 - j2 ) )
       and
            (        fixp    ( j3 + j4 - j5 ) and j3 + j4 - j5 >= 0
              or not numberp ( j3 + j4 - j5 ) )
       and
            (        fixp    ( j4 + j5 - j3 ) and j4 + j5 - j3 >= 0
              or not numberp ( j4 + j5 - j3 ) )
       and
            (        fixp    ( j5 + j3 - j4 ) and j5 + j3 - j4 >= 0
              or not numberp ( j5 + j3 - j4 ) )
       and
            (        fixp    ( j1 + j3 - j6 ) and j1 + j3 - j6 >= 0
              or not numberp ( j1 + j3 - j6 ) )
       and
            (        fixp    ( j3 + j6 - j1 ) and j3 + j6 - j1 >= 0
              or not numberp ( j3 + j6 - j1 ) )
       and
            (        fixp    ( j6 + j1 - j3 ) and j6 + j1 - j3 >= 0
              or not numberp ( j6 + j1 - j3 ) )
       and
            (        fixp    ( j2 + j4 - j6 ) and j2 + j4 - j6 >= 0
              or not numberp ( j2 + j4 - j6 ) )
       and
            (        fixp    ( j4 + j6 - j2 ) and j4 + j6 - j2 >= 0
              or not numberp ( j4 + j6 - j2 ) )
       and
            (        fixp    ( j6 + j2 - j4 ) and j6 + j2 - j4 >= 0
              or not numberp ( j6 + j2 - j4 ) )
     )
    then flag takes affirmed
    else flag takes denied $
    triangular!_cond!_check!_mode!_ takes Racah $
      Racah!_trang!_cond!_1!_save!_ takes j1 + j2 - j5 $
      Racah!_trang!_cond!_2!_save!_ takes j2 + j5 - j1 $
      Racah!_trang!_cond!_3!_save!_ takes j5 + j1 - j2 $
      Racah!_trang!_cond!_4!_save!_ takes j3 + j4 - j5 $
      Racah!_trang!_cond!_5!_save!_ takes j4 + j5 - j3 $
      Racah!_trang!_cond!_6!_save!_ takes j5 + j3 - j4 $
      Racah!_trang!_cond!_7!_save!_ takes j1 + j3 - j6 $
      Racah!_trang!_cond!_8!_save!_ takes j3 + j6 - j1 $
      Racah!_trang!_cond!_9!_save!_ takes j6 + j1 - j3 $
      Racah!_trang!_cond!_a!_save!_ takes j2 + j4 - j6 $
      Racah!_trang!_cond!_b!_save!_ takes j4 + j6 - j2 $
      Racah!_trang!_cond!_c!_save!_ takes j6 + j2 - j4 $
      Racah!_trang!_cond!_d!_save!_ takes j1 $
      Racah!_trang!_cond!_e!_save!_ takes j2 $
      Racah!_trang!_cond!_f!_save!_ takes j3 $
      Racah!_trang!_cond!_g!_save!_ takes j4 $
      Racah!_trang!_cond!_h!_save!_ takes j5 $
      Racah!_trang!_cond!_i!_save!_ takes j6 $
    return flag $
    end $
%
%
% Calculate a given Racah coefficient to give its concrete form.
%
procedure ww!#0!_ ( index!_, j1, j2, j3, j4, j5, j6 ) $
  begin scalar wrk0, wrk1, wrk2, wrk3, wrk4, wrk5 $
    if    not numberp j1 or not numberp j2 or not numberp j3
       or not numberp j4 or not numberp j5 or not numberp j6
       then do
      begin $
        clear p!_, q!_ $
        reduction!_control!_ ( 1 ) $
        wrk0 takes waux!#00!_ ( j1, j2, j3, j4, j5, j6 ) $
        on exp $
        wrk1 takes wrk0 $
        wrk0 takes waux!#01!_ ( j1, j2, j3, j4, j5, j6 ) $
        on exp $
        wrk2 takes wrk0 $
        wrk3 takes wrk1 * wrk2 $
        wrk0 takes wrk2 * wrk3 $
        reduction!_control!_ ( 2 ) $
        on exp $
        wrk1 takes wrk0 $
        reduction!_control!_ ( 3 ) $
        on exp $
        wrk4 takes sqrt num wrk1 $
        wrk5 takes sqrt den wrk1 $
        wrk3 takes hatch!_half!_eggs!_ ( wrk5 ) $
        wrk2 takes hatch!_half!_eggs!_ ( wrk4 ) $
        wrk0 takes wrk2 / wrk3 $
        if index!_ is 0 then wrk2 takes                  wrk0 $
        if index!_ is 1 then wrk2 takes 
          ( - 1 ) ** ( - p!_ * ( j2 + j3 - j5 - j6 ) ) * wrk0 $
        if index!_ is 2 then wrk2 takes 
          ( - 1 ) ** ( - p!_ * ( j1 + j4 - j5 - j6 ) ) * wrk0 $
        end
       else do
      begin $
        p!_ takes 1 $
        wrk0 takes sqrt ( waux!#00!_ ( j1, j2, j3, j4, j5, j6 ) )
                        * waux!#01!_ ( j1, j2, j3, j4, j5, j6 )   $
        if index!_ is 0 then wrk2 takes wrk0 $
        if index!_ is 1 then wrk2 takes 
          ( - 1 ) ** ( - p!_ * ( j2 + j3 - j5 - j6 ) ) * wrk0 $
        if index!_ is 2 then wrk2 takes 
          ( - 1 ) ** ( - p!_ * ( j1 + j4 - j5 - j6 ) ) * wrk0 $
        clear p!_ $
        if phase!_ is abstract and wrk2 isnot 0 then
           if part ( dummy1!_ * dummy2!_ * wrk2 , 0 ) is minus then
              wrk2 takes wrk2 / ( - 1 ) * ( - 1 ) ** ( - p!_ ) $
        end $
    return wrk2 $
    end $
%
procedure ww!#0!_2 ( index!_, j1, j2, j3, j4, j5, j6 ) $
  begin scalar wrk0, wrk1, wrk2, wrk3 $
    if    not numberp j1 or not numberp j2 or not numberp j3
       or not numberp j4 or not numberp j5 or not numberp j6
       then do
      begin $
        clear p!_, q!_ $
        reduction!_control!_ ( 1 ) $
        wrk0 takes waux!#00!_ ( j1, j2, j3, j4, j5, j6 ) $
        on exp $
        wrk1 takes wrk0 $
        wrk0 takes waux!#01!_ ( j1, j2, j3, j4, j5, j6 ) $
        on exp $
        wrk2 takes wrk0 $
        wrk3 takes wrk1 * wrk2 $
        wrk0 takes wrk2 * wrk3 $
        reduction!_control!_ ( 2 ) $
        on exp $
        wrk1 takes wrk0 $
        reduction!_control!_ ( 3 ) $
        on exp $
        wrk0 takes wrk1 $
        if index!_ is 0 then wrk1 takes                    wrk0 $
        if index!_ is 1 then wrk1 takes 
          ( - 1 )**( - p!_ * ( j2 + j3 - j5 - j6 ) * 2 ) * wrk0 $
        if index!_ is 2 then wrk1 takes 
          ( - 1 )**( - p!_ * ( j1 + j4 - j5 - j6 ) * 2 ) * wrk0 $
        wrk3 takes sub ( p!_ = q!_, wrk1 ) $
        wrk2 takes sub ( q!_ = p!_, wrk3 ) $
        end
       else do
      begin $
        p!_ takes 1 $
        wrk0 takes        waux!#00!_ ( j1, j2, j3, j4, j5, j6 )
                        * waux!#01!_ ( j1, j2, j3, j4, j5, j6 ) ** 2  $
        if index!_ is 0 then wrk2 takes                    wrk0 $
        if index!_ is 1 then wrk2 takes 
          ( - 1 )**( - p!_ * ( j2 + j3 - j5 - j6 ) * 2 ) * wrk0 $
        if index!_ is 2 then wrk2 takes 
          ( - 1 )**( - p!_ * ( j1 + j4 - j5 - j6 ) * 2 ) * wrk0 $
        clear p!_ $
        if phase!_ is abstract and wrk2 isnot 0 then
           if part ( dummy1!_ * dummy2!_ * wrk2 , 0 ) is minus then
              wrk2 takes wrk2 / ( - 1 ) * ( -1 ) ** ( - p!_ ) $
        end $
    return wrk2 $
    end $
%
%
% Calculate a Racah coefficient when one of the arguments is zero.
%
procedure ww!#1!_ ( index!_, j1, j2, j3, j4, j5, j6 ) $
    begin scalar wrk0, wrk1, wrk2 $
      wrk0 takes   ( - 1 ) ** ( p!_ * ( j6 - j2 - j4 ) )
              * krdlta ( j1, j2 )
              * krdlta ( j3, j4 )
              / sqrt ( egg ( 2 * j2 + 1 ) )
              / sqrt ( egg ( 2 * j4 + 1 ) ) $
      if index!_ is 0 then wrk1 takes                  wrk0 $
      if index!_ is 1 then wrk1 takes 
        ( - 1 ) ** ( - p!_ * ( j2 + j3 - j5 - j6 ) ) * wrk0 $
      if index!_ is 2 then wrk1 takes 
        ( - 1 ) ** ( - p!_ * ( j1 + j4 - j5 - j6 ) ) * wrk0 $
      return wrk1 $
      end $
%
procedure ww!#1!_2 ( index!_, j1, j2, j3, j4, j5, j6 ) $
    begin scalar wrk0, wrk1, wrk2 $
      wrk0 takes   ( - 1 ) ** ( p!_ * ( ( j6 - j2 - j4 ) * 2 ) )
              * krdlta ( j1, j2 ) ** 2
              * krdlta ( j3, j4 ) ** 2
              / egg ( 2 * j2 + 1 )
              / egg ( 2 * j4 + 1 ) $
      if index!_ is 0 then wrk1 takes                    wrk0 $
      if index!_ is 1 then wrk1 takes 
        ( - 1 )**( - p!_ * ( j2 + j3 - j5 - j6 ) * 2 ) * wrk0 $
      if index!_ is 2 then wrk1 takes 
        ( - 1 )**( - p!_ * ( j1 + j4 - j5 - j6 ) * 2 ) * wrk0 $
      return wrk1 $
      end $
%
%
operator w   $
operator w!_ $
%
operator w2   $
operator w2!_ $
%
operator w6j   $
operator w6j!_ $
%
operator w6j2   $
operator w6j2!_ $
%
%
% Racah coefficients.
%
% For general cases.
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 1
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#0!_ ( 0, j1, j2, j3, j4, j5, j6 ) $
%
% For special cases that one of j1 through j6 is zero.
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 11
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_ ( 2, j6, j3, j2, j5, j1, j4 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 12
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_ ( 1, j4, j6, j5, j1, j2, j3 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 13
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_ ( 1, j4, j5, j6, j1, j3, j2 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 14
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_ ( 2, j6, j2, j3, j5, j4, j1 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 15
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_ ( 0, j1, j2, j3, j4, j5, j6 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 16
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_ ( 0, j4, j2, j3, j1, j6, j5 ) $
%
% For cases that violate the triangular condition.
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 0
  let w!_ ( j1, j2, j3, j4, j5, j6 ) be 0 $
%
%
% Square of Racah coefficients
%
% For general cases.
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 1
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#0!_2 ( 0, j1, j2, j3, j4, j5, j6 ) $
%
% For special cases that one of j1 through j6 is zero.
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 11
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_2 ( 2, j6, j3, j2, j5, j1, j4 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 12
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_2 ( 1, j4, j6, j5, j1, j2, j3 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 13
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_2 ( 1, j4, j5, j6, j1, j3, j2 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 14
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_2 ( 2, j6, j2, j3, j5, j4, j1 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 15
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_2 ( 0, j1, j2, j3, j4, j5, j6 ) $
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 16
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 
    ww!#1!_2 ( 0, j4, j2, j3, j1, j6, j5 ) $
%
% For cases that violate the triangular condition.
%
for all j1, j2, j3, j4, j5, j6 such that
    rcoefficient!_index!_ is 0
  let w2!_ ( j1, j2, j3, j4, j5, j6 ) be 0 $
%
%
% Find out the most optimum condition for calculating a Racah
% coefficient, and perform the calculation.
%
procedure racah!_aux!_proc!_ 
  ( racah!_aux!_index!_, j1, j2, j3, j4, j5, j6 ) $
    begin scalar wrk0, wrk1, trcsave, flag $
      racah!_reduction!_condition!_flag!_ takes affirmed $
      trcsave takes triangular!_cond!_check!_mode!_ $
      rcoefficient!_indexw!_ takes undetermined $
      rcoefficient!_index!_  takes undetermined $
      r1!_ takes j1 $ r2!_ takes j2 $ r3!_ takes j3 $
      r4!_ takes j4 $ r5!_ takes j5 $ r6!_ takes j6 $
      wrk0 takes 1 $
      flag takes tragcondrc ( j1, j2, j3, j4, j5, j6 ) $
      if flag is affirmed then do
        begin $
          if r1!_ is 0 then rcoefficient!_index!_ takes 11  $
          if r2!_ is 0 then rcoefficient!_index!_ takes 12  $
          if r3!_ is 0 then rcoefficient!_index!_ takes 13  $
          if r4!_ is 0 then rcoefficient!_index!_ takes 14  $
          if r5!_ is 0 then rcoefficient!_index!_ takes 15  $
          if r6!_ is 0 then rcoefficient!_index!_ takes 16  $
          if rcoefficient!_index!_ is undetermined then do
            begin $
              ww!_x!_ takes 1 $
              r1!_ takes make!_window!_for!_x!_ ( j1 ) $
              r2!_ takes make!_window!_for!_x!_ ( j2 ) $
              r3!_ takes make!_window!_for!_x!_ ( j3 ) $
              r4!_ takes make!_window!_for!_x!_ ( j4 ) $
              r5!_ takes make!_window!_for!_x!_ ( j5 ) $
              r6!_ takes make!_window!_for!_x!_ ( j6 ) $
              wrk0 takes ww!_x!_ $
              clear kpdel!_m!_, kpmin!_m!_, kpmax!_m!_ $
              rcoefficient!_index!_ takes 
                racah!_gatein!_check!_aux!_ 
                (  1, r1!_, r2!_, r3!_, r4!_, r5!_, r6!_ ) $
              kpdel!_ takes kpdel!_m!_ $ 
              kpmin!_ takes kpmin!_m!_ $ 
              kpmax!_ takes kpmax!_m!_ $
              end $
          end
                                                           else do
        begin $
          rcoefficient!_index!_ takes 0 $
          end $
    if rcoefficient!_index!_ isnot undetermined then do
      begin $
        on exp $
        reduction!_control!_ ( 0 ) $
        if racah!_aux!_index!_ is 1 then
          wrk2 takes wrk0 *  w!_    
            ( r1!_, r2!_, r3!_, r4!_, r5!_, r6!_ ) $
        if racah!_aux!_index!_ is 2 then
          wrk2 takes wrk0 *  w2!_   
            ( r1!_, r2!_, r3!_, r4!_, r5!_, r6!_ ) $
        if racah!_aux!_index!_ is 3 then
          wrk2 takes wrk0 *  w6j!_  
            ( r1!_, r2!_, r3!_, r4!_, r5!_, r6!_ ) $
        if racah!_aux!_index!_ is 4 then
          wrk2 takes wrk0 *  w6j2!_ 
            ( r1!_, r2!_, r3!_, r4!_, r5!_, r6!_ ) $
        reduction!_control!_ ( 8 ) $
        end
                                                else do
      begin $
        wrk1 takes program!_error $
        end $
    off exp $
    if racah!_call!_from!_x!_flag!_ is affirmed
      then wrk1 takes           wrk2 
      else wrk1 takes hatcheggs wrk2 $
    off exp $
    triangular!_cond!_check!_mode!_ takes trcsave $
    return wrk1 $
    end $
%
procedure racah!_gatein!_check!_aux!_ 
  ( index!_no!_, j1, j2, j3, j4, j5, j6 ) $
  begin scalar flag $
    flag takes findsumrangeracah ( j1, j2, j3, j4, j5, j6 ) $
    if flag is successful then do
      begin $
        if    not fixp kpdel!_m!_ and fixp kpdel!_
           or fixp kpdel!_m!_ and kpdel!_m!_ > kpdel!_ then do
          begin $
            kpdel!_m!_ takes kpdel!_ $
            kpmin!_m!_ takes kpmin!_ $
            kpmax!_m!_ takes kpmax!_ $
            rcoefficient!_indexw!_ takes index!_no!_ $
            end $
      end $
    return rcoefficient!_indexw!_ $
    end $
%
%
%  End of Auxiliary functions.
%
%
%
% 4.2. Racah coefficients.
%
%
for all j1, j2, j3, j4, j5, j6 such that
          symbolic!_mode!_ is inactive
    and   rcc!_red!_mode!_00 is active
    and ( j1 is 0 or j2 is 0 or j3 is 0 
       or j4 is 0 or j5 is 0 or j6 is 0
       or fixp ( J4 + J3 + J2 + J1 )
       or fixp ( J6 + J5 + J4 + J1 )
       or fixp ( J6 + J5 + J3 + J2 )
       or fixp ( - J5 + J4 + J3 ) or fixp ( J5 - J4 + J3 ) 
       or fixp ( J5 + J4 - J3 )
       or fixp ( - J5 + J2 + J1 ) or fixp ( J5 - J2 + J1 ) 
       or fixp ( J5 + J2 - J1 )
       or fixp ( - J6 + J4 + J2 ) or fixp ( J6 - J4 + J2 ) 
       or fixp ( J6 + J4 - J2 )
       or fixp ( - J6 + J3 + J1 ) or fixp ( J6 - J3 + J1 ) 
       or fixp ( J6 + J3 - J1 ) )
  let w    ( j1, j2, j3, j4, j5, j6 ) be
    begin scalar wrk0, wrk1, wrk2 $
      clear p!_ $
      wrk0 takes racah!_aux!_proc!_ ( 1, j1, j2, j3, j4, j5, j6 ) $
      wrk1 takes wrk0 $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
%
% 4.3. Squre of Racah coefficients.
%
%
for all j1, j2, j3, j4, j5, j6 such that
          symbolic!_mode!_ is inactive
    and   rcc!_red!_mode!_00 is active
    and ( j1 is 0 or j2 is 0 or j3 is 0 
       or j4 is 0 or j5 is 0 or j6 is 0
       or fixp ( J4 + J3 + J2 + J1 )
       or fixp ( J6 + J5 + J4 + J1 )
       or fixp ( J6 + J5 + J3 + J2 )
       or fixp ( - J5 + J4 + J3 ) or fixp ( J5 - J4 + J3 ) 
       or fixp ( J5 + J4 - J3 )
       or fixp ( - J5 + J2 + J1 ) or fixp ( J5 - J2 + J1 ) 
       or fixp ( J5 + J2 - J1 )
       or fixp ( - J6 + J4 + J2 ) or fixp ( J6 - J4 + J2 ) 
       or fixp ( J6 + J4 - J2 )
       or fixp ( - J6 + J3 + J1 ) or fixp ( J6 - J3 + J1 ) 
       or fixp ( J6 + J3 - J1 ) )
  let w2   ( j1, j2, j3, j4, j5, j6 ) be
    begin scalar wrk0, wrk1, wrk2 $
      clear p!_ $
      wrk0 takes racah!_aux!_proc!_ ( 2, j1, j2, j3, j4, j5, j6 ) $
      wrk1 takes wrk0 $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
%
% 4.4. Wigner's 6-j symbols.
%
%
for all j1, j2, j3, j4, j5, j6 such that
          symbolic!_mode!_ is inactive
    and   rcc!_red!_mode!_00 is active
    and ( j1 is 0 or j2 is 0 or j5 is 0 
       or j4 is 0 or j3 is 0 or j6 is 0
       or fixp ( J4 + j5 + J2 + J1 )
       or fixp ( J6 + j3 + J4 + J1 )
       or fixp ( J6 + j3 + j5 + J2 )
       or fixp ( - j3 + J4 + j5 ) or fixp ( j3 - J4 + j5 ) 
       or fixp ( j3 + J4 - j5 )
       or fixp ( - j3 + J2 + J1 ) or fixp ( j3 - J2 + J1 ) 
       or fixp ( j3 + J2 - J1 )
       or fixp ( - J6 + J4 + J2 ) or fixp ( J6 - J4 + J2 ) 
       or fixp ( J6 + J4 - J2 )
       or fixp ( - J6 + j5 + J1 ) or fixp ( J6 - j5 + J1 ) 
       or fixp ( J6 + j5 - J1 ) )
  let w6j  ( j1, j2, j3, j4, j5, j6 ) be
    begin scalar wrk0, wrk1, wrk2 $
      clear p!_ $
      wrk0 takes racah!_aux!_proc!_ ( 1, j1, j2, j5, j4, j3, j6 ) $
      wrk1 takes   wrk0 * ( - 1 ) ** ( p!_ * ( j1 + j2 + j5 + j4 ) ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
%
% 4.5. Squre of Wigner's 6-j symbols.
%
%
for all j1, j2, j3, j4, j5, j6 such that
          symbolic!_mode!_ is inactive
    and   rcc!_red!_mode!_00 is active
    and ( j1 is 0 or j2 is 0 or j5 is 0 
       or j4 is 0 or j3 is 0 or j6 is 0
       or fixp ( J4 + j5 + J2 + J1 )
       or fixp ( J6 + j3 + J4 + J1 )
       or fixp ( J6 + j3 + j5 + J2 )
       or fixp ( - j3 + J4 + j5 ) or fixp ( j3 - J4 + j5 ) 
       or fixp ( j3 + J4 - j5 )
       or fixp ( - j3 + J2 + J1 ) or fixp ( j3 - J2 + J1 ) 
       or fixp ( j3 + J2 - J1 )
       or fixp ( - J6 + J4 + J2 ) or fixp ( J6 - J4 + J2 ) 
       or fixp ( J6 + J4 - J2 )
       or fixp ( - J6 + j5 + J1 ) or fixp ( J6 - j5 + J1 ) 
       or fixp ( J6 + j5 - J1 ) )
  let w6j2  ( j1, j2, j3, j4, j5, j6 ) be
    begin scalar wrk0, wrk1, wrk2 $
      clear p!_ $
      wrk0 takes racah!_aux!_proc!_ ( 2, j1, j2, j5, j4, j3, j6 ) $
      wrk1 takes 
             wrk0 * ( - 1 ) ** ( p!_ * ( j1 + j2 + j5 + j4 ) * 2 ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
%
%
% 4.6. Transformation between the variety of notations
%
% Notes: (0)   The following features are effective only when
%            the variable symbolic!_mode!_ is set active.
%            This mode can be set by calling a procedure evaluation, as
%
%              evaluation deactivate ;
%
%            The default is "activate", with which the reductions
%            will always be performed for the coefficients and symbols,
%            if possible.
%
%        (1)   When the variable w!_to!_w6j!_ is set "active"
%            and the variable w6j!_to!_w!_ is set "inactive"
%            at the call of w or w2, they are transformed into
%            the corresponding w6j or w6j2 prior to the
%            reduction of the coefficients. This mode can be
%            set by calling a procedure notation, as
%
%              notation Wigner ;
%
%            The default is "Asitis", with which no transformation
%            will be performed between the Racah's and 6j-symbols.
%
%        (2)   When the variable w!_to!_w6j!_ is set "inactive"
%            and the variable w6j!_to!_w!_ is set "active"
%            at the call of  w6j or w6j2, they are transformed into
%            the corresponding w or w2 prior to the
%            reduction of the coefficients. This mode can be
%            set by calling a procedure notation, as
%
%              notation Rose ;
%
%            The default is "Asitis", with which NO transformation
%            will be performed between the Racah's and 6j-symbols.
%
%        (3)   When both the variables  w!_to!_w6j!_ and w6j!_to!_w!_
%            are set "inactive", NO transformation will be performed.
%            This mode is the default and can be reset by calling a
%            procedure notation, as
%
%              notation Asitis ;
%
%        (4)   When both the variables  w!_to!_w6j!_ and w6j!_to!_w!_
%            are set "active", NO transformation will be performed
%            at all.
%
%
for all j1, j2, j3, j4, j5, j6 such that w!_to!_w6j!_ is   active
                                     and w6j!_to!_w!_ is inactive
                                     and symbolic!_mode!_ is active
  let w ( j1, j2, j3, j4, j5, j6 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                       w6j ( j1, j2, j5, j4, j3, j6 )
                     * ( - 1 ) ** ( p!_ * ( j1 + j2 + j3 + j4 ) ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
for all j1, j2, j3, j4, j5, j6 such that w!_to!_w6j!_ is inactive
                                     and w6j!_to!_w!_ is   active
                                     and symbolic!_mode!_ is active
  let w6j ( j1, j2, j5, j4, j3, j6 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                 w ( j1, j2, j3, j4, j5, j6 )
               / ( - 1 ) ** ( p!_ * ( j1 + j2 + j3 + j4 ) ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
for all j1, j2, j3, j4, j5, j6 such that w!_to!_w6j!_ is   active
                                     and w6j!_to!_w!_ is inactive
                                     and symbolic!_mode!_ is active
  let w2 ( j1, j2, j3, j4, j5, j6 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                 w6j2 ( j1, j2, j5, j4, j3, j6 )
               * ( - 1 ) ** ( p!_ * ( ( j1 + j2 + j3 + j4 ) * 2 ) ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
for all j1, j2, j3, j4, j5, j6 such that w!_to!_w6j!_ is inactive
                                     and w6j!_to!_w!_ is   active
                                     and symbolic!_mode!_ is active
  let w6j2( j1, j2, j5, j4, j3, j6 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      wrk1 takes
                 w2 ( j1, j2, j3, j4, j5, j6 )
               / ( - 1 ) ** ( p!_ * ( ( j1 + j2 + j3 + j4 ) * 2 ) ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
%
%
% End of Chapter 4. Racah coefficients and Wigner's 6j-Symbols
%
%
%
%
% Chapter 5. 9j-symbols or X-coefficients and related coefficients
%
%
% The full path name of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter50.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter50.r,v 1.1 92/04/16 20:37:53 koikef Locked $'/
%
%
% 5.1. 9j-symbols or X-coefficients
%
%
% Auxiliary functions.
%
%
% Find out if the given set of arguments satisfies the triangular
% conditions. Returns "affirmed" when satisfied, and returns
% denied if unsatisfied.
%
% Stores triangular conditions into the memory.
%
procedure tragcond9j ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) $
  begin scalar flag $
    if
     (
            (        fixp    ( j1 + j2 - j3 ) and j1 + j2 - j3 >= 0
              or not numberp ( j1 + j2 - j3 ) )
       and
            (        fixp    ( j2 + j3 - j1 ) and j2 + j3 - j1 >= 0
              or not numberp ( j2 + j3 - j1 ) )
       and
            (        fixp    ( j3 + j1 - j2 ) and j3 + j1 - j2 >= 0
              or not numberp ( j3 + j1 - j2 ) )
       and
            (        fixp    ( j4 + j5 - j6 ) and j4 + j5 - j6 >= 0
              or not numberp ( j4 + j5 - j6 ) )
       and
            (        fixp    ( j5 + j6 - j4 ) and j5 + j6 - j4 >= 0
              or not numberp ( j5 + j6 - j4 ) )
       and
            (        fixp    ( j6 + j4 - j5 ) and j6 + j4 - j5 >= 0
              or not numberp ( j6 + j4 - j5 ) )
       and
            (        fixp    ( j7 + j8 - j9 ) and j7 + j8 - j9 >= 0
              or not numberp ( j7 + j8 - j9 ) )
       and
            (        fixp    ( j8 + j9 - j7 ) and j8 + j9 - j7 >= 0
              or not numberp ( j8 + j9 - j7 ) )
       and
            (        fixp    ( j9 + j7 - j8 ) and j9 + j7 - j8 >= 0
              or not numberp ( j9 + j7 - j8 ) )
       and
            (        fixp    ( j1 + j4 - j7 ) and j1 + j4 - j7 >= 0
              or not numberp ( j1 + j4 - j7 ) )
       and
            (        fixp    ( j4 + j7 - j1 ) and j4 + j7 - j1 >= 0
              or not numberp ( j4 + j7 - j1 ) )
       and
            (        fixp    ( j7 + j1 - j4 ) and j7 + j1 - j4 >= 0
              or not numberp ( j7 + j1 - j4 ) )
       and
            (        fixp    ( j2 + j5 - j8 ) and j2 + j5 - j8 >= 0
              or not numberp ( j2 + j5 - j8 ) )
       and
            (        fixp    ( j5 + j8 - j2 ) and j5 + j8 - j2 >= 0
              or not numberp ( j5 + j8 - j2 ) )
       and
            (        fixp    ( j8 + j2 - j5 ) and j8 + j2 - j5 >= 0
              or not numberp ( j8 + j2 - j5 ) )
       and
            (        fixp    ( j3 + j6 - j9 ) and j3 + j6 - j9 >= 0
              or not numberp ( j3 + j6 - j9 ) )
       and
            (        fixp    ( j6 + j9 - j3 ) and j6 + j9 - j3 >= 0
              or not numberp ( j6 + j9 - j3 ) )
       and
            (        fixp    ( j9 + j3 - j6 ) and j9 + j3 - j6 >= 0
              or not numberp ( j9 + j3 - j6 ) )
     )
    then do
      begin $
        triangular!_cond!_check!_mode!_ takes Xcoef $
          Xcoef!_trang!_cond!_1!_save!_ takes j1 + j2 - j3 $
          Xcoef!_trang!_cond!_2!_save!_ takes j2 + j3 - j1 $
          Xcoef!_trang!_cond!_3!_save!_ takes j3 + j1 - j2 $
          Xcoef!_trang!_cond!_4!_save!_ takes j4 + j5 - j6 $
          Xcoef!_trang!_cond!_5!_save!_ takes j5 + j6 - j4 $
          Xcoef!_trang!_cond!_6!_save!_ takes j6 + j4 - j5 $
          Xcoef!_trang!_cond!_7!_save!_ takes j7 + j8 - j9 $
          Xcoef!_trang!_cond!_8!_save!_ takes j8 + j9 - j7 $
          Xcoef!_trang!_cond!_9!_save!_ takes j9 + j7 - j8 $
          Xcoef!_trang!_cond!_a!_save!_ takes j1 + j4 - j7 $
          Xcoef!_trang!_cond!_b!_save!_ takes j4 + j7 - j1 $
          Xcoef!_trang!_cond!_c!_save!_ takes j7 + j1 - j4 $
          Xcoef!_trang!_cond!_d!_save!_ takes j2 + j5 - j8 $
          Xcoef!_trang!_cond!_e!_save!_ takes j5 + j8 - j2 $
          Xcoef!_trang!_cond!_f!_save!_ takes j8 + j2 - j5 $
          Xcoef!_trang!_cond!_g!_save!_ takes j3 + j6 - j9 $
          Xcoef!_trang!_cond!_h!_save!_ takes j6 + j9 - j3 $
          Xcoef!_trang!_cond!_i!_save!_ takes j9 + j3 - j6 $
          Xcoef!_trang!_cond!_j!_save!_ takes j1 $
          Xcoef!_trang!_cond!_k!_save!_ takes j2 $
          Xcoef!_trang!_cond!_l!_save!_ takes j3 $
          Xcoef!_trang!_cond!_m!_save!_ takes j4 $
          Xcoef!_trang!_cond!_n!_save!_ takes j5 $
          Xcoef!_trang!_cond!_o!_save!_ takes j6 $
          Xcoef!_trang!_cond!_p!_save!_ takes j7 $
          Xcoef!_trang!_cond!_q!_save!_ takes j8 $
          Xcoef!_trang!_cond!_r!_save!_ takes j9 $
          flag takes affirmed
        end
    else  flag takes denied $
    return flag $
    end $
%
%
% Find out the best combinations of lower and upper limits
% of the summation. Returns successful when found and returns
% unsuccessful when not be found.
%
procedure findsumrangexaux ( tm1!_, tm2!_ ) $
  begin scalar wrk1, flag $
    flag takes unsuccessful $
    if fixp ( tm2!_ - tm1!_ ) then do
      begin $
        wrk1 takes tm2!_ - tm1!_ $
        if wrk1 >= 0 then do
          begin $
            if fixp tdlx!_ and tdlx!_ > wrk1 or not fixp tdlx!_  
              then do
              begin $
                tmnx!_ takes    tm1!_ $
                tmxx!_ takes    tm2!_ $
                tdlx!_ takes tmxx!_ - tmnx!_ $
                end $
            end $
        flag takes successful $
        end $
    if fixp ( tm2!_ + tm1!_ ) then do
      begin $
        wrk1 takes tm2!_ + tm1!_ $
        if wrk1 >= 0 then do
          begin $
            if fixp tdlx!_ and tdlx!_ > wrk1 or not fixp tdlx!_  
              then do
              begin $
                tmnx!_ takes  - tm1!_ $
                tmxx!_ takes    tm2!_ $
                tdlx!_ takes tmxx!_ - tmnx!_ $
                end $
            end $
        flag takes successful $
        end $
    return flag $
    end $
%
%
% Find out the t-sum range for a given combination on the
% angular momentum values or variables.
%
% Returns successful if found, and returns unsuccessful
% if failed.
%
procedure findsumrangex ( j2, j3, j4, j6, j7, j8 ) $
  begin scalar flag $
    clear tdlx!_, tmnx!_, tmxx!_ $
    flag takes findsumrangexaux ( j3 - j7, j3 + j7 ) $
    flag takes findsumrangexaux ( j3 - j7, j2 + j4 ) $
    flag takes findsumrangexaux ( j3 - j7, j6 + j8 ) $
    flag takes findsumrangexaux ( j2 - j4, j3 + j7 ) $
    flag takes findsumrangexaux ( j2 - j4, j2 + j4 ) $
    flag takes findsumrangexaux ( j2 - j4, j6 + j8 ) $
    flag takes findsumrangexaux ( j6 - j8, j3 + j7 ) $
    flag takes findsumrangexaux ( j6 - j8, j2 + j4 ) $
    flag takes findsumrangexaux ( j6 - j8, j6 + j8 ) $
    if fixp tdlx!_ then findsumrangex takes   successful
                   else findsumrangex takes unsuccessful $
    return findsumrangex $
    end $
%
%
% Main body of the x-coefficient calculations.
%
% Perform the t-sum of
%
%
% (2t+1)w(j2,j4,j3,j7;t,j1)w(j4,j2,j6,j8;t,j5)w(j7,j3,j8,j6;t,j9),
%
% and multiply an appropriate phase according to the symmetry
% condition. 
%
procedure xw!#0!_ ( index!_, j1, j2, j3, j4, j5, j6, j7, j8, j9 ) $
  begin scalar wrk0, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, 
               opsave, phsave $
    racah!_call!_from!_x!_flag!_ takes affirmed $
    if    not numberp j1 or not numberp j2 or not numberp j3
       or not numberp j4 or not numberp j5 or not numberp j6
       or not numberp j7 or not numberp j8 or not numberp j9
    then do
      begin $
        clear p!_ $
        opsave takes out!_put!_form!_flag!_ $
        phsave takes phase!_ $
        outputform intermediate $
        phasemode concrete $
        tx!_ takes 0 $
        wrk0 takes 0 $
        while tx!_ <= tdlx!_ do
          begin $
            wrk4 takes egg ( 2 * ( tx!_ + tmnx!_ ) + 1 ) $
            if message!_level!_ is chatty then
              write "Processing W ( ",j2,", ", j4,", ", j3,", ",
              j7,", ", tx!_ + tmnx!_,", ", j1," ) ** 2."            $
            racah!_reduction!_condition!_flag!_ takes denied $
            on exp $
            wrk1 takes w2 ( j2, j4, j3, j7, tx!_ + tmnx!_, j1 ) $
            if racah!_reduction!_condition!_flag!_ is denied then
              wrk1 takes w ( j2, j4, j3, j7, tx!_ + tmnx!_, j1 ) ** 2 $
            on exp $
            clear p!_ $
            if     wrk1 isnot 0 
               and part ( dummy1!_ * dummy2!_ * wrk1, 0 ) is minus
              then wrk1 takes wrk1 / ( - 1 ) * ( - 1 ) ** p!_ $
            if message!_level!_ is chatty then
              write "The result is ", wrk1, " ." $
            if message!_level!_ is chatty then
             write "Processing W ( ",j4,", ", j2,", ", j6,", ",
             j8,", ", tx!_ + tmnx!_,", ", j5," ) ** 2."            $
            racah!_reduction!_condition!_flag!_ takes denied $
            on exp $
            wrk2 takes w2 ( j4, j2, j6, j8, tx!_ + tmnx!_, j5 ) $
            if racah!_reduction!_condition!_flag!_ is denied then
              wrk2 takes w ( j4, j2, j6, j8, tx!_ + tmnx!_, j5 ) ** 2 $
            on exp $
            clear p!_ $
            if     wrk2 isnot 0 
               and part ( dummy1!_ * dummy2!_ * wrk2, 0 ) is minus
              then wrk2 takes wrk2 / ( - 1 ) * ( - 1 ) ** p!_ $
            if message!_level!_ is chatty then
              write "The result is ", wrk2, " ." $
            if message!_level!_ is chatty then
             write "Processing W ( ",j7,", ", j3,", ", j8,", ",
             j6,", ", tx!_ + tmnx!_,", ", j9," ) ** 2."            $
            racah!_reduction!_condition!_flag!_ takes denied $
            on exp $
            wrk3 takes w2 ( j7, j3, j8, j6, tx!_ + tmnx!_, j9 ) $
            if racah!_reduction!_condition!_flag!_ is denied then
              wrk3 takes w ( j7, j3, j8, j6, tx!_ + tmnx!_, j9 ) ** 2 $
            on exp $
            clear p!_ $
            if     wrk3 isnot 0 
               and part ( dummy1!_ * dummy2!_ * wrk3, 0 ) is minus
              then wrk3 takes wrk3 / ( - 1 ) * ( - 1 ) ** p!_ $
            if message!_level!_ is chatty then
              write "The result is ", wrk3, " ." $
            on exp $
            wrk5 takes wrk1 * wrk2 * wrk3 $
            wrk2 takes num wrk5 $
            wrk3 takes den wrk5 $
            wrk5 takes sqrt wrk2 $
            wrk6 takes sqrt wrk3 $
            wrk2 takes xeval!_aux!_ ( wrk5, wrk6 ) $
            p!_ takes 1 $
            wrk0 takes wrk0 + wrk2 * wrk4 $
            tx!_ takes tx!_ + 1 $
            end $
        clear p!_ $
        if     wrk0 isnot 0 
           and part ( dummy1!_ * dummy2!_ * wrk0, 0 ) is minus
          then wrk0 takes wrk0 / ( - 1 ) * ( - 1 ) ** p!_ $
        phasemode phsave $
        wrk1 takes
             ( - 1 ) ** ( p!_ * (
             ( j1 + j2 + j3 + j4 + j5 + j6 + j7 + j8 + j9 ) 
             * ( 1 - index!_ ) )
             ) * wrk0 $
        if phase!_ is concrete then p!_ takes 1 $
        outputform opsave $
        reduction!_control!_ ( 0 ) $
        reduction!_control!_ ( 2 ) $
        reduction!_control!_ ( 3 ) $
        wrk2 takes wrk1 $
        reduction!_control!_ ( 8 ) $
        clear p!_ $
        end
    else do
      begin $
        phsave takes phase!_ $
        phasemode concrete $
        on exp $
        tx!_ takes 0 $
        wrk0 takes 0 $
        while tx!_ <= tdlx!_ do
          begin $
            wrk3 takes 2 * ( tx!_ + tmnx!_ ) + 1 $
            wrk4 takes w ( j2, j4, j3, j7, tx!_ + tmnx!_, j1 ) $
            wrk5 takes w ( j4, j2, j6, j8, tx!_ + tmnx!_, j5 ) $
            wrk6 takes w ( j7, j3, j8, j6, tx!_ + tmnx!_, j9 ) $
            wrk0 takes wrk0 + wrk3 * wrk4 * wrk5 * wrk6 $
            tx!_ takes tx!_ + 1 $
            end $
        on exp $
        wrk2 takes
             ( - 1 ) ** ( 1 * (
             ( j1 + j2 + j3 + j4 + j5 + j6 + j7 + j8 + j9 ) 
             * ( 1 - index!_ ) )
             ) * wrk0 $
        clear p!_ $
        phasemode phsave $
        if phase!_ is abstract then do
          begin $
            if     wrk2 isnot 0 
               and part ( dummy1!_ * dummy!_ * wrk2, 0 ) is minus then
               wrk2 takes wrk2 * ( - 1 ) * ( - 1 ) ** ( p!_ ) $
            end $
        end $
    racah!_call!_from!_x!_flag!_ takes denied $
    off exp $
    wrk1 takes hatcheggs ( wrk2 ) $
    off exp $
    return wrk1 $
    end $
%
%
% An auxiliary function to the previously described
% procedure xw!#0!_. Makes the expressions neat.
% 
procedure xeval!_aux!_ ( x, y ) $
  begin scalar wrk1, wrk2, wrk3 $
    on exp $
    wrk2 takes hatch!_half!_eggs!_ ( x ) $
    wrk3 takes hatch!_half!_eggs!_ ( y ) $
    wrk1 takes wrk2 / wrk3 $
    work!_for!_meta!_mode!_of!_krdlta!_ takes wrk1 $
    meta!_mode!_of!_krdlta!_ takes active $
    wrk2 takes wrk1 $
    meta!_mode!_of!_krdlta!_ takes inactive $
    wrk2 takes work!_for!_meta!_mode!_of!_krdlta!_ $
    return wrk2 $
    end $
%
%
% The group of gateway routines that connect the 
% user interface routines and the procedure xw!#0!_.
% These operators are prepared to avoid tedious 
% permutations of the input angular momentum variables
% in the procedure xw!#0!_, and also to return 0
% immediately when the reduction conditions are 
% not fulfilled.
%
operator x   $
operator x!_ $
%
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
        xcoefficient!_index!_ is 1
  let x!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    xw!#0!_ ( 0, j1, j2, j3, j4, j5, j6, j7, j8, j9 ) $
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
        xcoefficient!_index!_ is 2
  let x!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    xw!#0!_ ( 1, j2, j1, j3, j5, j4, j6, j8, j7, j9 ) $
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
        xcoefficient!_index!_ is 3
  let x!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    xw!#0!_ ( 1, j1, j3, j2, j4, j6, j5, j7, j9, j8 ) $
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
        xcoefficient!_index!_ is 4
  let x!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    xw!#0!_ ( 1, j3, j2, j1, j6, j5, j4, j9, j8, j7 ) $
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
        xcoefficient!_index!_ is 5
  let x!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    xw!#0!_ ( 0, j2, j3, j1, j5, j6, j4, j8, j9, j7 ) $
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
        xcoefficient!_index!_ is 6
  let x!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    xw!#0!_ ( 0, j3, j1, j2, j6, j4, j5, j9, j7, j8 ) $
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
        xcoefficient!_index!_ is 0
  let x!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be 0 $
%
%
% Find one of the best conditions for the calculation of a given
% 9j-symbol, and perform the calculation.
%
% The calculation is so controlled as to minimize the range of
% the t-summation in the procedure xw!#0!_.
%
% If the condition is determined, the control will be transferred to the
% group of gateway routines x!_  and returns the result of calculation.
% The reduction condition must be determined in this phase of the 
% calculation. If turnd out to be undetermined, we would have hit 
% a bug of the program.
%
procedure xcoef!_aux!_proc!_ ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) $
    begin scalar wrk0, wrk1, flag $
      xcoefficient!_indexw!_ takes undetermined $
      xcoefficient!_index!_  takes undetermined $
      j1!_ takes j1 $ j2!_ takes j2 $ j3!_ takes j3 $
      j4!_ takes j4 $ j5!_ takes j5 $ j6!_ takes j6 $
      j7!_ takes j7 $ j8!_ takes j8 $ j9!_ takes j9 $
      wrk0 takes 1 $
      flag takes tragcond9j ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) $
      if flag is affirmed then do
        begin $
          ww!_x!_ takes 1 $
          j1!_ takes make!_window!_for!_x!_ ( j1 ) $
          j2!_ takes make!_window!_for!_x!_ ( j2 ) $
          j3!_ takes make!_window!_for!_x!_ ( j3 ) $
          j4!_ takes make!_window!_for!_x!_ ( j4 ) $
          j5!_ takes make!_window!_for!_x!_ ( j5 ) $
          j6!_ takes make!_window!_for!_x!_ ( j6 ) $
          j7!_ takes make!_window!_for!_x!_ ( j7 ) $
          j8!_ takes make!_window!_for!_x!_ ( j8 ) $
          j9!_ takes make!_window!_for!_x!_ ( j9 ) $
          wrk0 takes ww!_x!_ $
          clear tdlx!_m!_, tmnx!_m!_, tmxx!_m!_ $
          xcoefficient!_index!_ takes xcoef!_gatein!_check!_aux!_ 
            (  1, j2!_, j3!_, j4!_, j6!_, j7!_, j8!_ ) $
          xcoefficient!_index!_ takes xcoef!_gatein!_check!_aux!_ 
           (  2, j1!_, j3!_, j5!_, j6!_, j8!_, j7!_ ) $
          xcoefficient!_index!_ takes xcoef!_gatein!_check!_aux!_ 
           (  3, j3!_, j2!_, j4!_, j5!_, j7!_, j9!_ ) $
          xcoefficient!_index!_ takes xcoef!_gatein!_check!_aux!_ 
           (  4, j2!_, j1!_, j6!_, j4!_, j9!_, j8!_ ) $
          xcoefficient!_index!_ takes xcoef!_gatein!_check!_aux!_ 
           (  5, j3!_, j1!_, j5!_, j4!_, j8!_, j9!_ ) $
          xcoefficient!_index!_ takes xcoef!_gatein!_check!_aux!_ 
           (  6, j1!_, j2!_, j6!_, j5!_, j9!_, j7!_ ) $
          tdlx!_ takes tdlx!_m!_ $ 
          tmnx!_ takes tmnx!_m!_ $ 
          tmxx!_ takes tmxx!_m!_ $
          end
        else do
        begin $
          xcoefficient!_index!_ takes 0 $
          end $
    if xcoefficient!_index!_ isnot undetermined then
        wrk1 takes wrk0 *  
        x!_ ( j1!_, j2!_, j3!_, j4!_, j5!_, j6!_, j7!_, j8!_, j9!_ )
                                                else
        wrk1 takes program!_error $
    return wrk1 $
    end $
%
%
% An auxiliary routine to find out the minimum of the t-sum range.
%
procedure xcoef!_gatein!_check!_aux!_ 
  ( index!_no!_, j2, j3, j4, j6, j7, j8 ) $
  begin scalar flag $
    flag takes findsumrangex ( j2, j3, j4, j6, j7, j8 ) $
    if flag is successful then do
      begin $
        if not fixp tdlx!_m!_ and fixp tdlx!_
            or fixp tdlx!_m!_ and tdlx!_m!_ > tdlx!_ then do
          begin $
            tdlx!_m!_ takes tdlx!_ $
            tmnx!_m!_ takes tmnx!_ $
            tmxx!_m!_ takes tmxx!_ $
            xcoefficient!_indexw!_ takes index!_no!_ $
            end $
      end $
    return xcoefficient!_indexw!_ $
    end $
%
%
% End of auxiliary functions.
%
%
% The user interface of the x-coefficients.
% Transfer the control to the substantial part
% of the reduction system when the input 
% turns out to be reductive. Otherwise,
% do nothing.
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
       symbolic!_mode!_ is inactive
   and xct!_red!_mode!_00 is active
   and (
       findsumrangex ( j2, j3, j4, j6, j7, j8 ) is successful
    or findsumrangex ( j1, j3, j5, j6, j8, j7 ) is successful
    or findsumrangex ( j3, j2, j4, j5, j7, j9 ) is successful
    or findsumrangex ( j2, j1, j6, j4, j9, j8 ) is successful
    or findsumrangex ( j3, j1, j5, j4, j8, j9 ) is successful
    or findsumrangex ( j1, j2, j6, j5, j9, j7 ) is successful
       )
  let x ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    begin wrk0 $
      wrk0 takes xcoef!_aux!_proc!_ 
        ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) $
      return wrk0 $
      end $
%
%
%
% 5.2. Coefficient for transformation between jj and LS representations
%
operator tjjls $
%
for all j1, j2, j3, j4, j5, j6, j7, j8, j9 such that
          symbolic!_mode!_ is inactive
  let tjjls ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) be
    begin scalar wrk1, wrk2 $
      clear p!_ $
      off exp $
      on gcd $ off ezgcd $
      wrk1 takes sqrt (   ( 2 * j3 + 1 ) * ( 2 * j6 + 1 )
                        * ( 2 * j7 + 1 ) * ( 2 * j8 + 1 ) )
                 * x ( j1, j2, j3, j4, j5, j6, j7, j8, j9 ) $
      if phase!_ is concrete then p!_ takes 1 $
      wrk2 takes wrk1 $
      clear p!_ $
      return wrk2 $
      end $
%
%
% End of Chapter 5. 9j-symbols or X-coefficients 
% and related coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% Appendix. Diagnostics routines.                                      %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% The full path name of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter98.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter98.r,v 1.1 92/04/16 20:37:55 koikef Locked $'/
%
%
% Appendix A1. Diagnostics routines
%
% A1.1 Unitarity check of the Rose's c-coefficients.
%
% Auxiliary function :
%
procedure cuofc!_aux ( j1, j2, m ) $
begin scalar j, m1,countx,county,jmin,jmax,m1min,m1max,mwdth,cwrk$
  for all n1, n2 such that numberp n1 and numberp n2
    let sqrt ( n1 ) * sqrt ( n2 ) be sqrt ( n1 * n2 ) $
  if m + j2 >  j1 then m1max takes     j1
                  else m1max takes m + j2 $
  if m - j2 < -j1 then m1min takes    -j1
                  else m1min takes m - j2 $
  mwdth takes m1max - m1min + 1 $
  matrix u(mwdth,mwdth),umx(mwdth,mwdth) $
    countx takes 0 $
    if abs( j2 - j1 ) > abs ( m ) then jmin takes abs ( j2 - j1 )
                                  else jmin takes abs ( m       ) $
    jmax takes j2 + j1 $
    j takes jmin $
    while j <= jmax do
      begin $
        countx takes countx + 1 $
        county takes 0 $
        m1 takes m1min $
        while m1 <= m1max do
          begin $
            county takes county + 1 $
            cwrk takes c ( j1, j2, j, m1, m-m1, m ) $
            u  ( countx, county ) takes cwrk $
            if countx  is  county then umx( countx, county ) takes 1 $
            if countx neq county then umx( countx, county ) takes 0 $
            m1 takes m1 + 1 $
          end $
        j takes j + 1 $
      end $
    on gcd $  off ezgcd $
    off exp$
    if ( u * tp ( u ) )  is  umx then
      write 
        "Case ",j1,", ",j2,", ",m," : ","OK! Thats unitary"
                                else
      write 
        "Case ",j1,", ",j2,", ",m," : ","Sorry! Thats NOT unitary" $
end $
%
%  Unitarity check of the Rose's c-coefficients.
%
procedure check!_unitarity!_of!_c ( j1max, j2max ) $
begin scalar j121, j221, j1wk, j2wk $
  j121 takes 2 * j1max + 1 $
  j221 takes 2 * j2max + 1 $
  if j1max >= 0 and j2max >= 0 then do
    begin $
      for j1wk takes 1 : j121 do
        for j2wk takes 1 : j221 do
          begin $
            j1 takes ( j1wk - 1 ) / 2 $
            j2 takes ( j2wk - 1 ) / 2 $
            mmin takes - ( j1 + j2 ) $
            mmax takes     j1 + j2   $
            m takes mmin $
            while m <= mmax do
            begin$
              cuofc!_aux ( j1, j2, m ) $
              m takes m + 1 $
            end$
          end $
    end
                               else do
    begin $
      write "Improper argument values." $
    end $
end $
%
%
%
% A1.2 Symmetry property check of the Rose's c-coefficients
%
% Auxiliary function :
%
procedure csofc!_aux ( j1, j2, m ) $
begin scalar j, m1, m2, jmin,jmax,m1min,m1max $
  for all n1, n2 such that numberp n1 and numberp n2
    let sqrt ( n1 ) * sqrt ( n2 ) is sqrt ( n1 * n2 ) $
  if m + j2 >  j1 then m1max takes     j1
                  else m1max takes m + j2 $
  if m - j2 < -j1 then m1min takes    -j1
                  else m1min takes m - j2 $
    if abs( j2 - j1 ) > abs ( m ) then jmin takes abs ( j2 - j1 )
                                  else jmin takes abs ( m       ) $
    jmax takes j2 + j1 $
    j takes jmin $
    while j <= jmax do
      begin $
        m1 takes m1min $
        while m1 <= m1max do
          begin $
            m2 takes m - m1 $
            cwrk takes  c ( j1, j2, j , m1, m2, m  )$
            if (     cwrk
                 neq
                     (-1) ** ( j1 + j2 - j )
                   * c ( j2, j1, j , m2, m1, m  )
                 or
                     cwrk
                 neq
                     (-1) ** ( j1 + j2 - j )
                   * c ( j1, j2, j ,-m1,-m2,-m  )
                 or
                     cwrk
                   neq
                     (-1) ** ( j1 - m1 ) * sqrt(2*j +1)/sqrt(2*j2+1)
                   * c ( j1, j , j2, m1,-m ,-m2 )
                 or
                     cwrk
                   neq
                     (-1) ** ( j2 + m2 ) * sqrt(2*j +1)/sqrt(2*j1+1)
                   * c ( j , j2, j1,-m , m2,-m1 )
                 or
                     cwrk
                   neq
                     (-1) ** ( j1 - m1 ) * sqrt(2*j +1)/sqrt(2*j2+1)
                   * c ( j , j1, j2, m ,-m1, m2 )
                 or
                     cwrk
                   neq
                     (-1) ** ( j2 + m2 ) * sqrt(2*j +1)/sqrt(2*j1+1)
                   * c ( j2, j , j1,-m2, m , m1 )
                )
            then
              write "Case ",j1,", ",j2,", ",j,", ",m1,", ",m,", : ",
                    "Sorry! Thats NOT symmetric."
            $
            m1 takes m1 + 1 $
          end $
          write "Case ",j1,", ",j2,", ",j,", ",m," ",
                    "has been processed." $

        j takes j + 1 $
      end $
  for all n1, n2 such that numberp n1 and numberp n2
    clear sqrt ( n1 ) * sqrt ( n2 ) $
end $
%
% Check symmetry of the Rose's c-coefficients :
%
procedure check!_symmetry!_of!_c ( j1max, j2max ) $
begin scalar j121, j221, j1wk, j2wk $
  j121 takes 2 * j1max + 1 $
  j221 takes 2 * j2max + 1 $
  if j1max >= 0 and j2max >= 0 then do
    begin $
      for j1wk takes 1 : j121 do
        for j2wk takes 1 : j221 do
          begin $
            j1 takes ( j1wk - 1 ) / 2 $
            j2 takes ( j2wk - 1 ) / 2 $
            mmin takes - ( j1 + j2 ) $
            mmax takes     j1 + j2   $
            m takes mmin $
            while m <= mmax do
            begin$
              csofc!_aux ( j1, j2, m ) $
              m takes m + 1 $
            end$
          end $
    end
                               else do
    begin $
      write "Improper argument values." $
    end $
end $
%
%
% A1.3 Unitarity check of transformation coefficients trpopa.
%
procedure check!_unitarity!_of!_trpopa ( n, m ) $
begin scalar count,cl,qmin,qmax,nwdth,lmin $
  phase!_ takes concrete $
  nwdth takes n - abs ( m ) $
  lmin takes abs ( m ) $
  matrix u(nwdth,nwdth),ut(nwdth,nwdth),
         s(nwdth,nwdth),umx(nwdth,nwdth) $
  qmin takes - ( n - abs(m) - 1 ) $
  qmax takes   ( n - abs(m) - 1 ) $
  cl takes 0 $
  for l takes lmin : n-1 do
    begin $
      cl    takes cl + 1 $
      count takes 0 $
      for q takes qmin step 2 until qmax do
        begin $
          count takes count + 1 $
          u ( cl, count ) takes trpopa ( n, q, l, m ) $
          ut( count, cl ) takes u ( cl, count ) $
          if cl  is  count then umx(cl,count) takes 1 $
          if cl neq count then umx(cl,count) takes 0 $
        end $
    end $
    s takes u * ut $
    if s is umx then write "OK! Thats unitary" $
end $
%
%
% End of Appendix. Diagnostics routines.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% Epilogue. Initialize global variables.                               %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% The full path name of the RCS file:
%
%  /'$Source: /private/mnt/users/koikef/analg/RCS/Chapter99.r,v $'/
%
% The version number of the current program:
%
%  /'$Header: Chapter99.r,v 1.1 92/04/16 20:37:58 koikef Locked $'/
%
%
% Initialize global variables for controlling the
% program execution.
%
% Note: Do NOT edit any lines that include symbols 
%       with under bars: "_". Any symbols with under bars
%       are for the system use only.
%
% Constants :
%
clear   active $
clear inactive $
%
clear   activate $
clear deactivate $
%
clear Rose $
clear Wigner $
clear Asitis $
%
clear abstract $
clear concrete $
%
clear reduce $
clear retain $
%
clear final $
clear intermediate $
%
clear show $
clear nop $
%
clear affirmed $
clear denied $
%
clear   successful $
clear unsuccessful $
%
clear fullinteger $
clear halfinteger $
clear typeconflict $
%
clear quiet  $
clear chatty $
%
clear p!_ $
clear q!_ $
clear zero $
%
clear dummy!_  $
clear dummy1!_ $
clear dummy2!_ $
%
clear sw0!_0, sw0!_1, sw0!_2, sw0!_3, sw0!_4, sw0!_5, sw0!_6, sw0!_7 $
%
clear sw1!_0, sw1!_1, sw1!_2, sw1!_3, sw1!_4, sw1!_5, sw1!_6, sw1!_7 $
%
clear sw2!_0, sw2!_1, sw2!_2, sw2!_3, sw2!_4, sw2!_5, sw2!_6, sw2!_7 $
%
clear sw3!_0, sw3!_1, sw3!_2, sw3!_3, sw3!_4, sw3!_5, sw3!_6, sw3!_7 $
%
%
% Run time switches :
%
defaultswitchposition 0 $
%
defaultswitchposition 1 $
%
defaultswitchposition 2 $
%
defaultswitchposition 3 $
%
%
% Runtime modes :
%
evaluation activate $
%
notation Asitis $
%
autocontrol activate $
%
indirectreduction deactivate $
%
phasemode concrete $
%
factred activate $
%
factdiv activate $
%
pfactred activate $
%
outputform final $
%
polynomialformforphfunc activate $
%
reduction!_control!_ 5 $
%
sineformforphfunc deactivate $
%
reduction!_control!_ 6 $
%
runtimereporting quiet $
%
%
% Runtime flags :
%
meta!_mode!_of!_krdlta!_ takes inactive $
%
racah!_call!_from!_x!_flag!_ takes denied $
%
%
% Runtime system flags :
%
off exp $
%
off ezgcd $
%
on mcd $
%
% End of Epilogue. Initialize global variables.
%
$ END $
