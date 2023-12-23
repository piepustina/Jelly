function grad_int_r_i_X_dr_i = grad_int_r_X_dr(in1,in2)
%grad_int_r_X_dr
%    grad_int_r_i_X_dr_i = grad_int_r_X_dr(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:50:06

L_0 = in2(1,:);
R = in2(2,:);
deltaL = in1(3,:);
phi = in1(2,:);
rho = in2(3,:);
theta = in1(1,:);
t2 = L_0+deltaL;
t3 = R.^2;
t4 = phi.^2;
t5 = theta.^2;
t6 = 1.0./L_0;
t7 = t2.^2;
t8 = t4+t5;
t9 = 1.0./t8;
t11 = sqrt(t8);
t10 = t9.^2;
t12 = 1.0./t11;
t16 = cos(t11);
t17 = sin(t11);
t18 = t11.*4.6910077030668e-2;
t19 = t11.*7.692346550528415e-1;
t20 = t11.*2.307653449471585e-1;
t21 = t11.*9.53089922969332e-1;
t32 = -t11;
t35 = t11./2.0;
t13 = t12.^3;
t14 = t12.^5;
t15 = t12.^7;
t22 = t16.^2;
t23 = t17.^2;
t24 = cos(t18);
t25 = cos(t19);
t26 = sin(t18);
t27 = sin(t19);
t28 = cos(t20);
t29 = cos(t21);
t30 = sin(t20);
t31 = sin(t21);
t33 = t16-1.0;
t34 = L_0.*t16.*2.0;
t36 = cos(t35);
t37 = sin(t35);
t56 = t12.*t17;
t70 = L_0.*t17.*t32;
t38 = -t34;
t39 = L_0.*t22.*2.0;
t40 = L_0.*t23.*2.0;
t41 = t24-1.0;
t42 = t25-1.0;
t43 = t28-1.0;
t44 = t29-1.0;
t45 = t36-1.0;
t46 = L_0.*t4.*t24.*4.6910077030668e-2;
t47 = L_0.*t5.*t24.*4.6910077030668e-2;
t48 = L_0.*t4.*t25.*7.692346550528415e-1;
t49 = L_0.*t5.*t25.*7.692346550528415e-1;
t50 = L_0.*t4.*t28.*2.307653449471585e-1;
t51 = L_0.*t5.*t28.*2.307653449471585e-1;
t52 = L_0.*t4.*t29.*9.53089922969332e-1;
t53 = L_0.*t5.*t29.*9.53089922969332e-1;
t54 = (L_0.*t4.*t36)./2.0;
t55 = (L_0.*t5.*t36)./2.0;
t62 = t12.*t26;
t63 = t12.*t27;
t64 = t12.*t30;
t65 = t12.*t31;
t72 = t24.*t32;
t73 = t25.*t32;
t74 = t28.*t32;
t75 = t29.*t32;
t76 = t2.*t56;
t83 = L_0.*t26.*t32;
t84 = L_0.*t27.*t32;
t85 = L_0.*t30.*t32;
t86 = L_0.*t31.*t32;
t87 = t12.*t37;
t92 = t32.*t36;
t93 = L_0.*t32.*t37;
t98 = phi.*t2.*t10.*t33.*theta;
t105 = t2.*t4.*t10.*t33;
t106 = t2.*t5.*t10.*t33;
t108 = phi.*t2.*t13.*t26.*theta.*4.6910077030668e-2;
t109 = phi.*t2.*t13.*t27.*theta.*7.692346550528415e-1;
t114 = phi.*t2.*t13.*t30.*theta.*2.307653449471585e-1;
t115 = phi.*t2.*t13.*t31.*theta.*9.53089922969332e-1;
t117 = t2.*t4.*t13.*t26.*4.6910077030668e-2;
t118 = t2.*t5.*t13.*t26.*4.6910077030668e-2;
t120 = t2.*t4.*t13.*t27.*7.692346550528415e-1;
t122 = t2.*t5.*t13.*t27.*7.692346550528415e-1;
t125 = t2.*t4.*t13.*t30.*2.307653449471585e-1;
t126 = t2.*t5.*t13.*t30.*2.307653449471585e-1;
t127 = t2.*t4.*t13.*t31.*9.53089922969332e-1;
t128 = t2.*t5.*t13.*t31.*9.53089922969332e-1;
t141 = (phi.*t2.*t13.*t37.*theta)./2.0;
t142 = t9.*t16.*t33;
t143 = (t2.*t4.*t13.*t37)./2.0;
t144 = (t2.*t5.*t13.*t37)./2.0;
t151 = -t2.*t13.*(t11-t17);
t155 = phi.*t2.*t14.*theta.*(t11-t17).*3.0;
t156 = t2.*t4.*t14.*(t11-t17).*3.0;
t157 = t2.*t5.*t14.*(t11-t17).*3.0;
t88 = t2.*t62;
t89 = t2.*t63;
t90 = t2.*t64;
t91 = t2.*t65;
t94 = t2.*t9.*t41;
t95 = t2.*t9.*t42;
t96 = t2.*t9.*t43;
t97 = t2.*t9.*t44;
t100 = t2.*t87;
t107 = t2.*t9.*t45;
t119 = phi.*t2.*t10.*t41.*theta.*2.0;
t121 = phi.*t2.*t10.*t42.*theta.*2.0;
t123 = phi.*t2.*t10.*t43.*theta.*2.0;
t124 = phi.*t2.*t10.*t44.*theta.*2.0;
t130 = t2.*t4.*t10.*t41.*2.0;
t131 = t2.*t4.*t10.*t42.*2.0;
t132 = t2.*t5.*t10.*t41.*2.0;
t133 = t2.*t5.*t10.*t42.*2.0;
t134 = t2.*t4.*t10.*t43.*2.0;
t135 = t2.*t4.*t10.*t44.*2.0;
t136 = t2.*t5.*t10.*t43.*2.0;
t137 = t2.*t5.*t10.*t44.*2.0;
t138 = phi.*t2.*t10.*t45.*theta.*2.0;
t139 = t2.*t4.*t10.*t45.*2.0;
t140 = t2.*t5.*t10.*t45.*2.0;
t145 = t17+t72;
t146 = t17+t73;
t147 = t17+t74;
t148 = t17+t75;
t149 = t2.*t142;
t150 = t17+t92;
t166 = t38+t39+t40+t50+t51+t70+t85;
t167 = t38+t39+t40+t52+t53+t70+t86;
t168 = t38+t39+t40+t46+t47+t70+t83;
t169 = t38+t39+t40+t48+t49+t70+t84;
t171 = t38+t39+t40+t54+t55+t70+t93;
t110 = -t94;
t111 = -t95;
t112 = -t96;
t113 = -t97;
t129 = -t107;
t162 = t98+t108+t119+t155;
t163 = t98+t109+t121+t155;
t164 = t98+t114+t123+t155;
t165 = t98+t115+t124+t155;
t170 = t98+t138+t141+t155;
t177 = phi.*t6.*t7.*t15.*t145.*t168.*theta;
t178 = phi.*t6.*t7.*t15.*t146.*t169.*theta;
t179 = phi.*t6.*t7.*t15.*t147.*t166.*theta;
t180 = phi.*t6.*t7.*t15.*t148.*t167.*theta;
t186 = phi.*t6.*t7.*t15.*t150.*t171.*theta;
t187 = t105+t110+t117+t130+t151+t156;
t188 = t105+t111+t120+t131+t151+t156;
t189 = t106+t110+t118+t132+t151+t157;
t190 = t106+t111+t122+t133+t151+t157;
t191 = t105+t112+t125+t134+t151+t156;
t192 = t105+t113+t127+t135+t151+t156;
t193 = t106+t112+t126+t136+t151+t157;
t194 = t106+t113+t128+t137+t151+t157;
t195 = t105+t129+t139+t143+t151+t156;
t196 = t106+t129+t140+t144+t151+t157;
t197 = -t165.*(-t76+t91-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17));
t198 = -t162.*(-t76+t88-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17));
t199 = -t163.*(-t76+t89-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17));
t200 = -t164.*(-t76+t90-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17));
t201 = -t170.*(-t76+t100-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17));
t202 = t177+t198;
t203 = t180+t197;
t204 = t178+t199;
t205 = t179+t200;
t210 = t186+t201;
t206 = rho.*t3.*t205.*pi.*4.786286704993665e-1;
t207 = rho.*t3.*t202.*pi.*2.369268850561891e-1;
t208 = rho.*t3.*t203.*pi.*2.369268850561891e-1;
t209 = rho.*t3.*t204.*pi.*4.786286704993665e-1;
t211 = rho.*t3.*t210.*pi.*5.688888888889778e-1;
t212 = t206+t207+t208+t209+t211;
t213 = (L_0.*t212)./2.0;
et1 = rho.*t3.*pi.*(t189.*(-t76+t88-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t5.*t6.*t7.*t15.*t145.*t168).*2.369268850561891e-1+rho.*t3.*pi.*(t190.*(-t76+t89-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t5.*t6.*t7.*t15.*t146.*t169).*4.786286704993665e-1+rho.*t3.*pi.*(t193.*(-t76+t90-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t5.*t6.*t7.*t15.*t147.*t166).*4.786286704993665e-1+rho.*t3.*pi.*(t194.*(-t76+t91-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t5.*t6.*t7.*t15.*t148.*t167).*2.369268850561891e-1;
et2 = rho.*t3.*pi.*(t196.*(-t76+t100-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t5.*t6.*t7.*t15.*t150.*t171).*5.688888888889778e-1;
et3 = rho.*t3.*pi.*(t187.*(-t76+t88-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t4.*t6.*t7.*t15.*t145.*t168).*2.369268850561891e-1+rho.*t3.*pi.*(t188.*(-t76+t89-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t4.*t6.*t7.*t15.*t146.*t169).*4.786286704993665e-1+rho.*t3.*pi.*(t191.*(-t76+t90-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t4.*t6.*t7.*t15.*t147.*t166).*4.786286704993665e-1+rho.*t3.*pi.*(t192.*(-t76+t91-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t4.*t6.*t7.*t15.*t148.*t167).*2.369268850561891e-1;
et4 = rho.*t3.*pi.*(t195.*(-t76+t100-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t4.*t6.*t7.*t15.*t150.*t171).*5.688888888889778e-1;
et5 = rho.*t3.*pi.*(phi.*t13.*t145.*(-t76+t88-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-phi.*t2.*t13.*t145.*(-t56+t62-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*2.369268850561891e-1+rho.*t3.*pi.*(phi.*t13.*t146.*(-t76+t89-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-phi.*t2.*t13.*t146.*(-t56+t63-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*4.786286704993665e-1;
et6 = rho.*t3.*pi.*(phi.*t13.*t147.*(-t76+t90-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-phi.*t2.*t13.*t147.*(-t56+t64-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*4.786286704993665e-1+rho.*t3.*pi.*(phi.*t13.*t148.*(-t76+t91-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-phi.*t2.*t13.*t148.*(-t56+t65-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*2.369268850561891e-1;
et7 = rho.*t3.*pi.*(phi.*t13.*t150.*(-t76+t100-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-phi.*t2.*t13.*t150.*(-t56+t87-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*5.688888888889778e-1;
et8 = rho.*t3.*pi.*(t13.*t145.*theta.*(-t76+t88-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t2.*t13.*t145.*theta.*(-t56+t62-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*2.369268850561891e-1+rho.*t3.*pi.*(t13.*t146.*theta.*(-t76+t89-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t2.*t13.*t146.*theta.*(-t56+t63-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*4.786286704993665e-1;
et9 = rho.*t3.*pi.*(t13.*t147.*theta.*(-t76+t90-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t2.*t13.*t147.*theta.*(-t56+t64-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*4.786286704993665e-1+rho.*t3.*pi.*(t13.*t148.*theta.*(-t76+t91-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t2.*t13.*t148.*theta.*(-t56+t65-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*2.369268850561891e-1;
et10 = rho.*t3.*pi.*(t13.*t150.*theta.*(-t76+t100-t149+t2.*t4.*t10.*t17.*(t11-t17)+t2.*t5.*t10.*t17.*(t11-t17))-t2.*t13.*t150.*theta.*(-t56+t87-t142+t4.*t10.*t17.*(t11-t17)+t5.*t10.*t17.*(t11-t17))).*5.688888888889778e-1;
mt1 = [-t213,(L_0.*(et3+et4))./2.0,(L_0.*(et5+et6+et7))./2.0,L_0.*(et1+et2).*(-1.0./2.0),t213,L_0.*(et8+et9+et10).*(-1.0./2.0),L_0.*(rho.*t3.*pi.*(phi.*t2.*t13.*t145.*t189-t2.*t13.*t145.*t162.*theta).*2.369268850561891e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t146.*t190-t2.*t13.*t146.*t163.*theta).*4.786286704993665e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t147.*t193-t2.*t13.*t147.*t164.*theta).*4.786286704993665e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t148.*t194-t2.*t13.*t148.*t165.*theta).*2.369268850561891e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t150.*t196-t2.*t13.*t150.*t170.*theta).*5.688888888889778e-1).*(-1.0./2.0)];
mt2 = [L_0.*(rho.*t3.*pi.*(phi.*t2.*t13.*t145.*t162-t2.*t13.*t145.*t187.*theta).*2.369268850561891e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t146.*t163-t2.*t13.*t146.*t188.*theta).*4.786286704993665e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t147.*t164-t2.*t13.*t147.*t191.*theta).*4.786286704993665e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t148.*t165-t2.*t13.*t148.*t192.*theta).*2.369268850561891e-1+rho.*t3.*pi.*(phi.*t2.*t13.*t150.*t170-t2.*t13.*t150.*t195.*theta).*5.688888888889778e-1).*(-1.0./2.0),0.0];
grad_int_r_i_X_dr_i = reshape([mt1,mt2],3,3);
end