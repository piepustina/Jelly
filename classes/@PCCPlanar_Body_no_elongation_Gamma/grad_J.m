function grad_J = grad_J(theta,in2)
%grad_J
%    grad_J = grad_J(THETA,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    13-Oct-2023 18:12:56

L_0 = in2(1,:);
R = in2(2,:);
rho = in2(3,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = theta.*3.125e-2;
t5 = theta.*6.25e-2;
t6 = theta.*1.25e-1;
t7 = theta.*2.5e-1;
t8 = theta.*5.545614591379636e-2;
t9 = theta.*1.109122918275927e-1;
t10 = theta.*2.218245836551854e-1;
t11 = theta.*4.436491673103708e-1;
t12 = theta.*7.043854086203644e-3;
t13 = theta.*1.408770817240729e-2;
t14 = theta.*2.817541634481458e-2;
t15 = theta.*5.635083268962916e-2;
t16 = R.^2;
t33 = 1.0./pi;
t43 = 1.0./rho;
t44 = theta./4.0;
t47 = theta./8.0;
t48 = theta./1.6e+1;
t49 = theta./3.2e+1;
t17 = cos(t4);
t18 = cos(t5);
t19 = cos(t6);
t20 = cos(t7);
t21 = sin(t4);
t22 = sin(t5);
t23 = sin(t6);
t24 = sin(t7);
t25 = cos(t8);
t26 = cos(t9);
t27 = cos(t10);
t28 = cos(t11);
t29 = cos(t12);
t30 = cos(t13);
t31 = cos(t14);
t32 = cos(t15);
t34 = 1.0./t16;
t35 = sin(t8);
t36 = sin(t9);
t37 = sin(t10);
t38 = sin(t11);
t39 = sin(t12);
t40 = sin(t13);
t41 = sin(t14);
t42 = sin(t15);
t45 = cos(t44);
t46 = sin(t44);
t50 = cos(t47);
t51 = cos(t48);
t52 = cos(t49);
t53 = L_0.*t3.*t27.*t28.*t36.*9.841229182759271e-2;
t54 = L_0.*t3.*t26.*t28.*t37.*1.968245836551854e-1;
t55 = L_0.*t3.*t26.*t27.*t38.*3.936491673103708e-1;
t56 = L_0.*t3.*t30.*t31.*t32.*-1.127016653792583e-1;
t57 = L_0.*t3.*t30.*t31.*t32.*1.127016653792583e-1;
t58 = L_0.*t3.*t31.*t32.*t40.*1.587708172407289e-3;
t59 = L_0.*t3.*t30.*t32.*t41.*3.175416344814578e-3;
t60 = L_0.*t3.*t30.*t31.*t42.*6.350832689629156e-3;
t61 = L_0.*t2.*t19.*t20.*t22.*3.125e-2;
t62 = L_0.*t2.*t18.*t20.*t23.*6.25e-2;
t63 = L_0.*t2.*t18.*t19.*t24.*1.25e-1;
t64 = L_0.*t3.*t18.*t19.*t20.*5.0e-1;
t65 = L_0.*t3.*t19.*t20.*t22.*3.125e-2;
t66 = L_0.*t3.*t18.*t20.*t23.*6.25e-2;
t67 = L_0.*t3.*t18.*t19.*t24.*1.25e-1;
t68 = L_0.*t2.*t26.*t27.*t28.*-8.872983346207417e-1;
t70 = L_0.*t2.*t27.*t28.*t36.*9.841229182759271e-2;
t71 = L_0.*t2.*t26.*t28.*t37.*1.968245836551854e-1;
t72 = L_0.*t2.*t26.*t27.*t38.*3.936491673103708e-1;
t73 = L_0.*t2.*t30.*t31.*t32.*-1.127016653792583e-1;
t75 = L_0.*t2.*t31.*t32.*t40.*1.587708172407289e-3;
t76 = L_0.*t2.*t30.*t32.*t41.*3.175416344814578e-3;
t77 = L_0.*t2.*t30.*t31.*t42.*6.350832689629156e-3;
t78 = L_0.*t3.*t26.*t27.*t28.*-8.872983346207417e-1;
t79 = L_0.*t3.*t26.*t27.*t28.*8.872983346207417e-1;
t80 = L_0.*rho.*t16.*t26.*t27.*t28.*1.548627738659433;
t81 = L_0.*rho.*t16.*t27.*t28.*t36.*1.717618516825001e-1;
t82 = L_0.*rho.*t16.*t26.*t28.*t37.*3.435237033650001e-1;
t83 = L_0.*rho.*t16.*t26.*t27.*t38.*6.870474067300002e-1;
t84 = L_0.*rho.*t16.*t30.*t31.*t32.*1.967015133348961e-1;
t85 = L_0.*rho.*t16.*t31.*t32.*t40.*2.771073516932898e-3;
t86 = L_0.*rho.*t16.*t30.*t32.*t41.*5.542147033865795e-3;
t87 = L_0.*rho.*t16.*t30.*t31.*t42.*1.108429406773159e-2;
t88 = L_0.*rho.*t16.*t18.*t19.*t20.*1.396263401595146;
t89 = L_0.*rho.*t16.*t19.*t20.*t22.*8.726646259969664e-2;
t90 = L_0.*rho.*t16.*t18.*t20.*t23.*1.745329251993933e-1;
t91 = L_0.*rho.*t16.*t18.*t19.*t24.*3.490658503987865e-1;
t92 = (L_0.*t2.*t18.*t19.*t20)./2.0;
t93 = L_0.*rho.*t16.*t30.*t31.*t32.*pi.*6.261203632181017e-2;
t94 = L_0.*rho.*t16.*t26.*t27.*t28.*pi.*4.929435192337454e-1;
t96 = L_0.*t2.*t25.*t26.*t27.*t28.*3.936491673103708e-1;
t97 = L_0.*t2.*t29.*t30.*t31.*t32.*6.350832689629156e-3;
t98 = L_0.*t3.*t25.*t26.*t27.*t28.*-3.936491673103708e-1;
t100 = L_0.*t2.*t17.*t18.*t19.*t20.*1.25e-1;
t101 = L_0.*t3.*t29.*t30.*t31.*t32.*-6.350832689629156e-3;
t103 = L_0.*t2.*t25.*t26.*t27.*t38.*-8.872983346207417e-1;
t104 = L_0.*t2.*t25.*t26.*t27.*t38.*8.872983346207417e-1;
t105 = L_0.*t2.*t26.*t27.*t35.*t38.*-4.920614591379636e-2;
t107 = L_0.*t2.*t25.*t27.*t36.*t38.*-9.841229182759271e-2;
t109 = L_0.*t2.*t25.*t26.*t37.*t38.*-1.968245836551854e-1;
t111 = L_0.*t3.*t17.*t18.*t19.*t20.*-1.25e-1;
t113 = L_0.*t2.*t29.*t30.*t31.*t42.*-1.127016653792583e-1;
t114 = L_0.*t2.*t29.*t30.*t31.*t42.*1.127016653792583e-1;
t115 = L_0.*t2.*t30.*t31.*t39.*t42.*-7.938540862036445e-4;
t117 = L_0.*t2.*t29.*t31.*t40.*t42.*-1.587708172407289e-3;
t119 = L_0.*t2.*t29.*t30.*t41.*t42.*-3.175416344814578e-3;
t121 = L_0.*t3.*t25.*t26.*t27.*t38.*-8.872983346207417e-1;
t123 = L_0.*t3.*t26.*t27.*t35.*t38.*4.920614591379636e-2;
t124 = L_0.*t3.*t25.*t27.*t36.*t38.*9.841229182759271e-2;
t125 = L_0.*t3.*t25.*t26.*t37.*t38.*1.968245836551854e-1;
t126 = L_0.*t2.*t18.*t19.*t21.*t24.*-1.5625e-2;
t128 = L_0.*t2.*t17.*t19.*t22.*t24.*-3.125e-2;
t130 = L_0.*t2.*t17.*t18.*t23.*t24.*-6.25e-2;
t132 = L_0.*t3.*t29.*t30.*t31.*t42.*-1.127016653792583e-1;
t134 = L_0.*t3.*t30.*t31.*t39.*t42.*7.938540862036445e-4;
t135 = L_0.*t3.*t29.*t31.*t40.*t42.*1.587708172407289e-3;
t136 = L_0.*t3.*t29.*t30.*t41.*t42.*3.175416344814578e-3;
t137 = L_0.*t3.*t18.*t19.*t21.*t24.*1.5625e-2;
t138 = L_0.*t3.*t17.*t19.*t22.*t24.*3.125e-2;
t139 = L_0.*t3.*t17.*t18.*t23.*t24.*6.25e-2;
t140 = L_0.*rho.*t16.*t17.*t18.*t19.*t20.*-3.490658503987865e-1;
t142 = L_0.*rho.*t16.*t17.*t18.*t19.*t24.*1.396263401595146;
t143 = L_0.*rho.*t16.*t18.*t19.*t21.*t24.*4.363323129984832e-2;
t144 = L_0.*rho.*t16.*t17.*t19.*t22.*t24.*8.726646259969664e-2;
t145 = L_0.*rho.*t16.*t17.*t18.*t23.*t24.*1.745329251993933e-1;
t146 = (L_0.*t2.*t17.*t18.*t19.*t24)./2.0;
t147 = (L_0.*t3.*t17.*t18.*t19.*t24)./2.0;
t148 = L_0.*rho.*t16.*t25.*t26.*t27.*t28.*-6.870474067300002e-1;
t150 = L_0.*rho.*t16.*t29.*t30.*t31.*t32.*-1.108429406773159e-2;
t152 = L_0.*rho.*t16.*t25.*t26.*t27.*t38.*1.548627738659433;
t153 = L_0.*rho.*t16.*t26.*t27.*t35.*t38.*8.588092584125003e-2;
t154 = L_0.*rho.*t16.*t25.*t27.*t36.*t38.*1.717618516825001e-1;
t155 = L_0.*rho.*t16.*t25.*t26.*t37.*t38.*3.435237033650001e-1;
t156 = L_0.*rho.*t16.*t29.*t30.*t31.*t42.*1.967015133348961e-1;
t157 = L_0.*rho.*t16.*t30.*t31.*t39.*t42.*1.385536758466449e-3;
t158 = L_0.*rho.*t16.*t29.*t31.*t40.*t42.*2.771073516932898e-3;
t159 = L_0.*rho.*t16.*t29.*t30.*t41.*t42.*5.542147033865795e-3;
t160 = L_0.*rho.*t16.*t29.*t30.*t31.*t42.*pi.*6.261203632181017e-2;
t161 = L_0.*rho.*t16.*t25.*t26.*t27.*t38.*pi.*4.929435192337454e-1;
t164 = (L_0.*t2.*t45.*t50.*t51)./2.0;
t165 = (L_0.*t3.*t45.*t50.*t51)./2.0;
t166 = L_0.*rho.*t16.*t45.*t50.*t51.*pi.*4.444444444443434e-1;
t168 = (L_0.*t2.*t46.*t50.*t51.*t52)./2.0;
t169 = (L_0.*t3.*t46.*t50.*t51.*t52)./2.0;
t170 = L_0.*rho.*t16.*t46.*t50.*t51.*t52.*pi.*4.444444444443434e-1;
t95 = -t92;
t162 = -t146;
t163 = -t147;
t171 = t80+t84+t88;
t172 = t93+t94+t166;
t173 = t142+t152+t156;
t179 = t160+t161+t170;
t187 = t81+t82+t83+t85+t86+t87+t89+t90+t91;
t198 = t140+t143+t144+t145+t148+t150+t153+t154+t155+t157+t158+t159;
t174 = t2.*t34.*t43.*t171.*1.591549430918953e-1;
t175 = t3.*t34.*t43.*t171.*-1.591549430918953e-1;
t177 = t2.*t34.*t43.*t173.*1.591549430918953e-1;
t178 = t3.*t34.*t43.*t173.*1.591549430918953e-1;
t180 = (t3.*t33.*t34.*t43.*t172)./2.0;
t181 = (t2.*t33.*t34.*t43.*t172)./2.0;
t183 = (t2.*t33.*t34.*t43.*t179)./2.0;
t184 = (t3.*t33.*t34.*t43.*t179)./2.0;
t188 = t2.*t34.*t43.*t187.*-1.591549430918953e-1;
t190 = t3.*t34.*t43.*t187.*-1.591549430918953e-1;
t199 = t2.*t34.*t43.*t198.*1.591549430918953e-1;
t200 = t3.*t34.*t43.*t198.*-1.591549430918953e-1;
t182 = -t181;
t185 = -t183;
t186 = -t184;
t202 = t53+t54+t55+t68+t96+t105+t107+t109+t121+t174+t178+t190+t199;
t203 = t70+t71+t72+t79+t98+t103+t123+t124+t125+t175+t177+t188+t200;
t204 = t58+t59+t60+t73+t97+t115+t117+t119+t132+t174+t178+t190+t199;
t205 = t57+t75+t76+t77+t101+t113+t134+t135+t136+t175+t177+t188+t200;
t206 = t61+t62+t63+t64+t111+t137+t138+t139+t162+t175+t177+t188+t200;
t207 = t65+t66+t67+t95+t100+t126+t128+t130+t163+t174+t178+t190+t199;
t192 = t78+t104+t180+t185;
t193 = t56+t114+t180+t185;
t197 = t164+t169+t182+t186;
t214 = rho.*t16.*t206.*pi.*(t165-t168-t180+t183).*8.888888888886868e-1;
t209 = rho.*t16.*t192.*t203.*pi.*5.555555555555556e-1;
t211 = rho.*t16.*t193.*t205.*pi.*5.555555555555556e-1;
t215 = -t214;
t216 = rho.*t16.*t197.*t207.*pi.*8.888888888886868e-1;
t218 = L_0.*(t209+t211+t215-t216+rho.*t16.*t202.*pi.*(t68+t121+t181+t184).*5.555555555555556e-1+rho.*t16.*t204.*pi.*(t73+t132+t181+t184).*5.555555555555556e-1).*(-1.0./2.0);
grad_J = reshape([(L_0.*(rho.*t16.*t203.*pi.*(t68+t121+t181+t184).*1.111111111111111+rho.*t16.*t205.*pi.*(t73+t132+t181+t184).*1.111111111111111-rho.*t16.*t197.*t206.*pi.*1.777777777777374))./2.0,0.0,t218,0.0,(L_0.*(rho.*t16.*pi.*(t203.*(t68+t121+t181+t184)+t192.*t202).*1.111111111111111+rho.*t16.*pi.*(t205.*(t73+t132+t181+t184)+t193.*t204).*1.111111111111111-rho.*t16.*pi.*(t197.*t206+t207.*(t165-t168-t180+t183)).*1.777777777777374))./2.0,0.0,t218,0.0,(L_0.*(rho.*t16.*t192.*t202.*pi.*1.111111111111111+rho.*t16.*t193.*t204.*pi.*1.111111111111111-rho.*t16.*t207.*pi.*(t165-t168-t180+t183).*1.777777777777374))./2.0],[3,3]);
end
