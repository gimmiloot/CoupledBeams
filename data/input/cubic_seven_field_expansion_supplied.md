# Кубическая семиполевая пространственная модель: полное разложение

Это результат аналитического раскрытия принятого в разговоре закона V^(0).
Файл не является изменением репозитория CoupledBeams, реализацией решателя
или физической валидацией нелинейной модели.

## Область и соглашения

Одно первоначально прямое однородное изотропное плечо, без внешних
распределённых сил и демпфирования. Коэффициенты могут различаться
между плечами, но постоянны внутри каждого плеча.
Сдвиговые жёсткости в обеих плоскостях одинаковы:
\(S_\parallel=S_\perp=S\).

Порядок полей в этом файле:
\[
Q=(u,w,v,\Phi,\psi,\theta,c).
\]
Материальная координата — \(s\), время — \(t\).
Все поля и их нормированные производные имеют первую амплитудную степень.
Геометрия и упругие коэффициенты при разложении фиксированы.

\[
a=(\Phi,-\psi,\theta)^T,\quad
R=\exp([a]_\times),\quad
\Gamma=R^T(e_1+U_s)-e_1,\quad U=(u,w,v)^T,
\]
\[
[\chi]_\times=R^TR_s,\qquad [\Omega]_\times=R^TR_t.
\]
\[
D_\Gamma=\operatorname{diag}(C,S,S),\quad
D_\chi=\operatorname{diag}(\mathcal C_T,B_\perp,B_\parallel),
\]
\[
J(c)=\operatorname{diag}\big(j_\perp+(1+c)^2j_\parallel,\ j_\perp,\ (1+c)^2j_\parallel\big).
\]
\[
\mathcal T=\frac m2|U_t|^2+\frac{j_\parallel}2c_t^2+
\frac12\Omega^TJ(c)\Omega,
\]
\[
\mathcal V^{(0)}=\frac12\Gamma^TD_\Gamma\Gamma+
\nu Cc\Gamma_1+\frac C2c^2+\frac H2c_s^2+\frac12\chi^TD_\chi\chi.
\]

Энергия V^(0) — принятое конститутивное допущение, не единственный
результат кинематики и не полная трёхмерная редукция.
Кинетическая энергия относится к принятому аффинному перемещению массы
сечения с масштабом 1+c, без дополнительной инерции депланации.

## Определение остатков

\[
\mathcal L_{\le4}=\mathcal T_2+\mathcal T_3+\mathcal T_4
-\mathcal V_2-\mathcal V_3-\mathcal V_4,
\]
\[
\mathcal E_{Q_j}=
\partial_t\frac{\partial\mathcal L_{\le4}}{\partial Q_{j,t}}+
\partial_s\frac{\partial\mathcal L_{\le4}}{\partial Q_{j,s}}-
\frac{\partial\mathcal L_{\le4}}{\partial Q_j}.
\]

Уравнение каждого поля:
\[
\mathcal E_{Q_j}^{[1]}+\mathcal E_{Q_j}^{[2]}+\mathcal E_{Q_j}^{[3]}=0.
\]
Это точная полиномиальная система действия L_(<=4); относительно исходной
неусечённой редуцированной модели отброшен остаток четвёртой и более
высокой амплитудной степени в уравнениях. Символы [1],[2],[3] означают
степень по полям, а не номер плеча. Никакие связи нерастяжимости,
зависимость поворотов от наклонов оси или c=-nu*u_s не наложены.

## Полностью раскрытые уравнения


### Поле $u$


\[
\begin{aligned}
\mathcal E_{u}^{[1]} &=m u_{tt} -C u_{ss} -C c_{s} \nu
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{u}^{[2]} &=C \psi \psi_{s} +C \theta \theta_{s} +S \psi v_{ss} \\
&\quad +S \psi_{s} v_{s} +S \theta w_{ss} +S \theta_{s} w_{s} \\
&\quad -C \psi v_{ss} -C \psi_{s} v_{s} -C \theta w_{ss} \\
&\quad -C \theta_{s} w_{s} -2 S \psi \psi_{s} -2 S \theta \theta_{s}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{u}^{[3]} &=C \psi^{2} u_{ss} +C \theta^{2} u_{ss} -S \psi^{2} u_{ss} \\
&\quad -S \theta^{2} u_{ss} +\frac{C \Phi \psi w_{ss}}{2} +\frac{C \Phi \psi_{s} w_{s}}{2} \\
&\quad +\frac{C \Phi_{s} \psi w_{s}}{2} +\frac{C c_{s} \nu \psi^{2}}{2} +\frac{C c_{s} \nu \theta^{2}}{2} \\
&\quad +\frac{\Phi S \theta v_{ss}}{2} +\frac{\Phi S \theta_{s} v_{s}}{2} +\frac{\Phi_{s} S \theta v_{s}}{2} \\
&\quad -2 S \psi \psi_{s} u_{s} -2 S \theta \theta_{s} u_{s} +2 C \psi \psi_{s} u_{s} \\
&\quad +2 C \theta \theta_{s} u_{s} -\frac{C \Phi \theta v_{ss}}{2} -\frac{C \Phi \theta_{s} v_{s}}{2} \\
&\quad -\frac{C \Phi_{s} \theta v_{s}}{2} -\frac{\Phi S \psi w_{ss}}{2} -\frac{\Phi S \psi_{s} w_{s}}{2} \\
&\quad -\frac{\Phi_{s} S \psi w_{s}}{2} +C c \nu \psi \psi_{s} +C c \nu \theta \theta_{s}
\end{aligned}
\]


### Поле $w$


\[
\begin{aligned}
\mathcal E_{w}^{[1]} &=S \theta_{s} +m w_{tt} -S w_{ss}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{w}^{[2]} &=S \theta u_{ss} +S \theta_{s} u_{s} -C \theta u_{ss} \\
&\quad -C \theta_{s} u_{s} -\frac{\Phi S \psi_{s}}{2} -\frac{\Phi_{s} S \psi}{2} \\
&\quad -C c \nu \theta_{s} -C c_{s} \nu \theta
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{w}^{[3]} &=S \theta^{2} w_{ss} +\frac{C \psi^{2} \theta_{s}}{2} -C \theta^{2} w_{ss} \\
&\quad -2 S \theta^{2} \theta_{s} -\frac{2 S \psi^{2} \theta_{s}}{3} -\frac{\Phi^{2} S \theta_{s}}{6} \\
&\quad +\frac{3 C \theta^{2} \theta_{s}}{2} +C \psi \psi_{s} \theta +S \psi \theta v_{ss} \\
&\quad +S \psi \theta_{s} v_{s} +S \psi_{s} \theta v_{s} +\frac{C \Phi \psi u_{ss}}{2} \\
&\quad +\frac{C \Phi \psi_{s} u_{s}}{2} +\frac{C \Phi_{s} \psi u_{s}}{2} -C \psi \theta v_{ss} \\
&\quad -C \psi \theta_{s} v_{s} -C \psi_{s} \theta v_{s} -2 C \theta \theta_{s} w_{s} \\
&\quad +2 S \theta \theta_{s} w_{s} -\frac{4 S \psi \psi_{s} \theta}{3} -\frac{\Phi S \psi u_{ss}}{2} \\
&\quad -\frac{\Phi S \psi_{s} u_{s}}{2} -\frac{\Phi_{s} S \psi u_{s}}{2} -\frac{\Phi \Phi_{s} S \theta}{3} \\
&\quad +\frac{C \Phi c \nu \psi_{s}}{2} +\frac{C \Phi c_{s} \nu \psi}{2} +\frac{C \Phi_{s} c \nu \psi}{2}
\end{aligned}
\]


### Поле $v$


\[
\begin{aligned}
\mathcal E_{v}^{[1]} &=S \psi_{s} +m v_{tt} -S v_{ss}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{v}^{[2]} &=S \psi u_{ss} +S \psi_{s} u_{s} +\frac{\Phi S \theta_{s}}{2} \\
&\quad +\frac{\Phi_{s} S \theta}{2} -C \psi u_{ss} -C \psi_{s} u_{s} \\
&\quad -C c \nu \psi_{s} -C c_{s} \nu \psi
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{v}^{[3]} &=S \psi^{2} v_{ss} +\frac{C \psi_{s} \theta^{2}}{2} -C \psi^{2} v_{ss} \\
&\quad -2 S \psi^{2} \psi_{s} -\frac{2 S \psi_{s} \theta^{2}}{3} -\frac{\Phi^{2} S \psi_{s}}{6} \\
&\quad +\frac{3 C \psi^{2} \psi_{s}}{2} +C \psi \theta \theta_{s} +S \psi \theta w_{ss} \\
&\quad +S \psi \theta_{s} w_{s} +S \psi_{s} \theta w_{s} +\frac{\Phi S \theta u_{ss}}{2} \\
&\quad +\frac{\Phi S \theta_{s} u_{s}}{2} +\frac{\Phi_{s} S \theta u_{s}}{2} -C \psi \theta w_{ss} \\
&\quad -C \psi \theta_{s} w_{s} -C \psi_{s} \theta w_{s} -2 C \psi \psi_{s} v_{s} \\
&\quad +2 S \psi \psi_{s} v_{s} -\frac{4 S \psi \theta \theta_{s}}{3} -\frac{C \Phi \theta u_{ss}}{2} \\
&\quad -\frac{C \Phi \theta_{s} u_{s}}{2} -\frac{C \Phi_{s} \theta u_{s}}{2} -\frac{\Phi \Phi_{s} S \psi}{3} \\
&\quad -\frac{C \Phi c \nu \theta_{s}}{2} -\frac{C \Phi c_{s} \nu \theta}{2} -\frac{C \Phi_{s} c \nu \theta}{2}
\end{aligned}
\]


### Поле $\Phi$


\[
\begin{aligned}
\mathcal E_{\Phi}^{[1]} &=\Phi_{tt} j_{\parallel} +\Phi_{tt} j_{\perp} -\mathcal C_T \Phi_{ss}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{\Phi}^{[2]} &=B_{\parallel} \psi_{s} \theta_{s} +j_{\perp} \psi_{t} \theta_{t} +\frac{B_{\parallel} \psi \theta_{ss}}{2} \\
&\quad +\frac{\mathcal C_T \psi_{ss} \theta}{2} +\frac{S \psi w_{s}}{2} +\frac{j_{\perp} \psi \theta_{tt}}{2} \\
&\quad -B_{\perp} \psi_{s} \theta_{s} -j_{\parallel} \psi_{t} \theta_{t} +2 \Phi_{t} c_{t} j_{\parallel} \\
&\quad +2 \Phi_{tt} c j_{\parallel} -\frac{B_{\perp} \psi_{ss} \theta}{2} -\frac{\mathcal C_T \psi \theta_{ss}}{2} \\
&\quad -\frac{S \theta v_{s}}{2} -\frac{j_{\parallel} \psi_{tt} \theta}{2}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{\Phi}^{[3]} &=\Phi_{tt} c^{2} j_{\parallel} +\frac{B_{\parallel} \Phi \psi_{s}^{2}}{2} +\frac{B_{\perp} \Phi \theta_{s}^{2}}{2} \\
&\quad -\frac{B_{\parallel} \Phi \theta_{s}^{2}}{2} -\frac{B_{\perp} \Phi \psi_{s}^{2}}{2} -\frac{\Phi j_{\parallel} \psi_{t}^{2}}{3} \\
&\quad -\frac{\Phi j_{\perp} \theta_{t}^{2}}{3} -\frac{\Phi_{tt} j_{\parallel} \theta^{2}}{3} -\frac{\Phi_{tt} j_{\perp} \psi^{2}}{3} \\
&\quad -\frac{B_{\parallel} \Phi_{ss} \psi^{2}}{4} -\frac{B_{\perp} \Phi_{ss} \theta^{2}}{4} -\frac{\mathcal C_T \Phi \psi_{s}^{2}}{6} \\
&\quad -\frac{\mathcal C_T \Phi \theta_{s}^{2}}{6} -\frac{\Phi S \psi^{2}}{12} -\frac{\Phi S \theta^{2}}{12} \\
&\quad -\frac{\Phi_{tt} j_{\parallel} \psi^{2}}{12} -\frac{\Phi_{tt} j_{\perp} \theta^{2}}{12} +\frac{\mathcal C_T \Phi_{ss} \psi^{2}}{3} \\
&\quad +\frac{\mathcal C_T \Phi_{ss} \theta^{2}}{3} +\frac{2 \Phi j_{\parallel} \theta_{t}^{2}}{3} +\frac{2 \Phi j_{\perp} \psi_{t}^{2}}{3} \\
&\quad +\frac{C \theta u_{s} v_{s}}{2} +\frac{S \psi u_{s} w_{s}}{2} -c j_{\parallel} \psi_{tt} \theta \\
&\quad -c_{t} j_{\parallel} \psi_{t} \theta -2 c j_{\parallel} \psi_{t} \theta_{t} +2 \Phi_{t} c c_{t} j_{\parallel} \\
&\quad -\frac{2 \Phi_{t} j_{\parallel} \theta \theta_{t}}{3} -\frac{2 \Phi_{t} j_{\perp} \psi \psi_{t}}{3} -\frac{B_{\parallel} \Phi_{s} \psi \psi_{s}}{2} \\
&\quad -\frac{B_{\perp} \Phi_{s} \theta \theta_{s}}{2} -\frac{C \psi u_{s} w_{s}}{2} -\frac{S \theta u_{s} v_{s}}{2} \\
&\quad -\frac{B_{\parallel} \Phi \theta \theta_{ss}}{6} -\frac{B_{\perp} \Phi \psi \psi_{ss}}{6} -\frac{\mathcal C_T \Phi \psi \psi_{ss}}{6} \\
&\quad -\frac{\mathcal C_T \Phi \theta \theta_{ss}}{6} -\frac{\Phi_{t} j_{\parallel} \psi \psi_{t}}{6} -\frac{\Phi_{t} j_{\perp} \theta \theta_{t}}{6} \\
&\quad -\frac{\Phi j_{\parallel} \psi \psi_{tt}}{12} -\frac{\Phi j_{\perp} \theta \theta_{tt}}{12} +\frac{\Phi S \psi v_{s}}{3} \\
&\quad +\frac{\Phi S \theta w_{s}}{3} +\frac{\Phi j_{\parallel} \theta \theta_{tt}}{3} +\frac{\Phi j_{\perp} \psi \psi_{tt}}{3} \\
&\quad +\frac{B_{\parallel} \Phi \psi \psi_{ss}}{4} +\frac{B_{\perp} \Phi \theta \theta_{ss}}{4} +\frac{2 \mathcal C_T \Phi_{s} \psi \psi_{s}}{3} \\
&\quad +\frac{2 \mathcal C_T \Phi_{s} \theta \theta_{s}}{3} +\frac{C c \nu \theta v_{s}}{2} -\frac{C c \nu \psi w_{s}}{2}
\end{aligned}
\]


### Поле $\psi$


\[
\begin{aligned}
\mathcal E_{\psi}^{[1]} &=S \psi +j_{\perp} \psi_{tt} -B_{\perp} \psi_{ss} \\
&\quad -S v_{s}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{\psi}^{[2]} &=C u_{s} v_{s} +\mathcal C_T \Phi_{s} \theta_{s} +\frac{B_{\perp} \Phi \theta_{ss}}{2} \\
&\quad +\frac{\mathcal C_T \Phi_{ss} \theta}{2} +\frac{\Phi S w_{s}}{2} +\frac{\Phi j_{\parallel} \theta_{tt}}{2} \\
&\quad -B_{\parallel} \Phi_{s} \theta_{s} -C \psi u_{s} -\Phi_{t} j_{\perp} \theta_{t} \\
&\quad -S u_{s} v_{s} +2 S \psi u_{s} -\frac{B_{\parallel} \Phi \theta_{ss}}{2} \\
&\quad -\frac{B_{\perp} \Phi_{ss} \theta}{2} -\frac{\Phi j_{\perp} \theta_{tt}}{2} -\frac{\Phi_{tt} j_{\parallel} \theta}{2} \\
&\quad +C c \nu v_{s} -C c \nu \psi
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{\psi}^{[3]} &=\frac{C \psi^{3}}{2} -\frac{2 S \psi^{3}}{3} +C \psi v_{s}^{2} \\
&\quad +S \psi u_{s}^{2} +\frac{B_{\parallel} \Phi_{s}^{2} \psi}{2} +\frac{C \psi \theta^{2}}{2} \\
&\quad +\frac{\mathcal C_T \psi \theta_{s}^{2}}{2} -C \psi u_{s}^{2} -S \psi v_{s}^{2} \\
&\quad +2 S \psi^{2} v_{s} -\frac{3 C \psi^{2} v_{s}}{2} -\frac{2 S \psi \theta^{2}}{3} \\
&\quad -\frac{B_{\parallel} \psi \theta_{s}^{2}}{2} -\frac{C \theta^{2} v_{s}}{2} -\frac{\mathcal C_T \Phi_{s}^{2} \psi}{2} \\
&\quad -\frac{j_{\perp} \psi \theta_{t}^{2}}{3} -\frac{\Phi^{2} j_{\perp} \psi_{tt}}{3} -\frac{B_{\parallel} \Phi^{2} \psi_{ss}}{4} \\
&\quad -\frac{\mathcal C_T \psi_{ss} \theta^{2}}{4} -\frac{B_{\perp} \Phi_{s}^{2} \psi}{6} -\frac{B_{\perp} \psi \theta_{s}^{2}}{6} \\
&\quad -\frac{\Phi^{2} S \psi}{12} -\frac{j_{\perp} \psi_{tt} \theta^{2}}{12} +\frac{B_{\perp} \Phi^{2} \psi_{ss}}{3} \\
&\quad +\frac{B_{\perp} \psi_{ss} \theta^{2}}{3} +\frac{\Phi^{2} j_{\parallel} \psi_{tt}}{4} +\frac{j_{\parallel} \psi_{tt} \theta^{2}}{4} \\
&\quad +\frac{\Phi^{2} S v_{s}}{6} +\frac{2 S \theta^{2} v_{s}}{3} +\frac{2 \Phi_{t}^{2} j_{\perp} \psi}{3} \\
&\quad +C \theta v_{s} w_{s} +\Phi c j_{\parallel} \theta_{tt} +\Phi c_{t} j_{\parallel} \theta_{t} \\
&\quad +\frac{\Phi \Phi_{t} j_{\parallel} \psi_{t}}{2} +\frac{\Phi S u_{s} w_{s}}{2} +\frac{j_{\parallel} \psi_{t} \theta \theta_{t}}{2} \\
&\quad -C \psi \theta w_{s} -\Phi_{t} c_{t} j_{\parallel} \theta -\Phi_{tt} c j_{\parallel} \theta \\
&\quad -S \theta v_{s} w_{s} -\frac{2 \Phi \Phi_{t} j_{\perp} \psi_{t}}{3} -\frac{B_{\parallel} \Phi \Phi_{s} \psi_{s}}{2} \\
&\quad -\frac{C \Phi u_{s} w_{s}}{2} -\frac{\mathcal C_T \psi_{s} \theta \theta_{s}}{2} -\frac{B_{\parallel} \psi \theta \theta_{ss}}{6} \\
&\quad -\frac{B_{\perp} \Phi \Phi_{ss} \psi}{6} -\frac{B_{\perp} \psi \theta \theta_{ss}}{6} -\frac{\mathcal C_T \Phi \Phi_{ss} \psi}{6} \\
&\quad -\frac{j_{\perp} \psi_{t} \theta \theta_{t}}{6} -\frac{\Phi \Phi_{tt} j_{\parallel} \psi}{12} -\frac{j_{\parallel} \psi \theta \theta_{tt}}{12} \\
&\quad -\frac{j_{\perp} \psi \theta \theta_{tt}}{12} +\frac{\Phi \Phi_{tt} j_{\perp} \psi}{3} +\frac{B_{\parallel} \Phi \Phi_{ss} \psi}{4} \\
&\quad +\frac{\mathcal C_T \psi \theta \theta_{ss}}{4} +\frac{2 B_{\perp} \Phi \Phi_{s} \psi_{s}}{3} +\frac{2 B_{\perp} \psi_{s} \theta \theta_{s}}{3} \\
&\quad +\frac{4 S \psi \theta w_{s}}{3} -C c \nu \psi u_{s} -\frac{C \Phi c \nu w_{s}}{2}
\end{aligned}
\]


### Поле $\theta$


\[
\begin{aligned}
\mathcal E_{\theta}^{[1]} &=S \theta +j_{\parallel} \theta_{tt} -B_{\parallel} \theta_{ss} \\
&\quad -S w_{s}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{\theta}^{[2]} &=B_{\perp} \Phi_{s} \psi_{s} +C u_{s} w_{s} +\Phi_{t} j_{\parallel} \psi_{t} \\
&\quad +\frac{B_{\parallel} \Phi_{ss} \psi}{2} +\frac{B_{\perp} \Phi \psi_{ss}}{2} +\frac{\Phi j_{\parallel} \psi_{tt}}{2} \\
&\quad +\frac{\Phi_{tt} j_{\perp} \psi}{2} -C \theta u_{s} -\mathcal C_T \Phi_{s} \psi_{s} \\
&\quad -S u_{s} w_{s} +2 S \theta u_{s} +2 c j_{\parallel} \theta_{tt} \\
&\quad +2 c_{t} j_{\parallel} \theta_{t} -\frac{B_{\parallel} \Phi \psi_{ss}}{2} -\frac{\mathcal C_T \Phi_{ss} \psi}{2} \\
&\quad -\frac{\Phi S v_{s}}{2} -\frac{\Phi j_{\perp} \psi_{tt}}{2} +C c \nu w_{s} \\
&\quad -C c \nu \theta
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{\theta}^{[3]} &=\frac{C \theta^{3}}{2} -\frac{2 S \theta^{3}}{3} +C \theta w_{s}^{2} \\
&\quad +S \theta u_{s}^{2} +c^{2} j_{\parallel} \theta_{tt} +\frac{B_{\perp} \Phi_{s}^{2} \theta}{2} \\
&\quad +\frac{C \psi^{2} \theta}{2} +\frac{\mathcal C_T \psi_{s}^{2} \theta}{2} -C \theta u_{s}^{2} \\
&\quad -S \theta w_{s}^{2} +2 S \theta^{2} w_{s} -\frac{3 C \theta^{2} w_{s}}{2} \\
&\quad -\frac{2 S \psi^{2} \theta}{3} -\frac{B_{\perp} \psi_{s}^{2} \theta}{2} -\frac{C \psi^{2} w_{s}}{2} \\
&\quad -\frac{\mathcal C_T \Phi_{s}^{2} \theta}{2} -\frac{j_{\parallel} \psi_{t}^{2} \theta}{3} -\frac{\Phi^{2} j_{\parallel} \theta_{tt}}{3} \\
&\quad -\frac{B_{\perp} \Phi^{2} \theta_{ss}}{4} -\frac{\mathcal C_T \psi^{2} \theta_{ss}}{4} -\frac{B_{\parallel} \Phi_{s}^{2} \theta}{6} \\
&\quad -\frac{B_{\parallel} \psi_{s}^{2} \theta}{6} -\frac{\Phi^{2} S \theta}{12} -\frac{j_{\parallel} \psi^{2} \theta_{tt}}{12} \\
&\quad +\frac{B_{\parallel} \Phi^{2} \theta_{ss}}{3} +\frac{B_{\parallel} \psi^{2} \theta_{ss}}{3} +\frac{\Phi^{2} j_{\perp} \theta_{tt}}{4} \\
&\quad +\frac{j_{\perp} \psi^{2} \theta_{tt}}{4} +\frac{\Phi^{2} S w_{s}}{6} +\frac{2 S \psi^{2} w_{s}}{3} \\
&\quad +\frac{2 \Phi_{t}^{2} j_{\parallel} \theta}{3} +C \psi v_{s} w_{s} +\Phi c j_{\parallel} \psi_{tt} \\
&\quad +\Phi c_{t} j_{\parallel} \psi_{t} +\frac{C \Phi u_{s} v_{s}}{2} +\frac{\Phi \Phi_{t} j_{\perp} \theta_{t}}{2} \\
&\quad +\frac{j_{\perp} \psi \psi_{t} \theta_{t}}{2} -C \psi \theta v_{s} -S \psi v_{s} w_{s} \\
&\quad +2 \Phi_{t} c j_{\parallel} \psi_{t} +2 c c_{t} j_{\parallel} \theta_{t} -\frac{2 \Phi \Phi_{t} j_{\parallel} \theta_{t}}{3} \\
&\quad -\frac{B_{\perp} \Phi \Phi_{s} \theta_{s}}{2} -\frac{\mathcal C_T \psi \psi_{s} \theta_{s}}{2} -\frac{\Phi S u_{s} v_{s}}{2} \\
&\quad -\frac{B_{\parallel} \Phi \Phi_{ss} \theta}{6} -\frac{B_{\parallel} \psi \psi_{ss} \theta}{6} -\frac{B_{\perp} \psi \psi_{ss} \theta}{6} \\
&\quad -\frac{\mathcal C_T \Phi \Phi_{ss} \theta}{6} -\frac{j_{\parallel} \psi \psi_{t} \theta_{t}}{6} -\frac{\Phi \Phi_{tt} j_{\perp} \theta}{12} \\
&\quad -\frac{j_{\parallel} \psi \psi_{tt} \theta}{12} -\frac{j_{\perp} \psi \psi_{tt} \theta}{12} +\frac{\Phi \Phi_{tt} j_{\parallel} \theta}{3} \\
&\quad +\frac{B_{\perp} \Phi \Phi_{ss} \theta}{4} +\frac{\mathcal C_T \psi \psi_{ss} \theta}{4} +\frac{2 B_{\parallel} \Phi \Phi_{s} \theta_{s}}{3} \\
&\quad +\frac{2 B_{\parallel} \psi \psi_{s} \theta_{s}}{3} +\frac{4 S \psi \theta v_{s}}{3} +\frac{C \Phi c \nu v_{s}}{2} \\
&\quad -C c \nu \theta u_{s}
\end{aligned}
\]


### Поле $c$


\[
\begin{aligned}
\mathcal E_{c}^{[1]} &=C c +c_{tt} j_{\parallel} -H c_{ss} \\
&\quad +C \nu u_{s}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{c}^{[2]} &=-\Phi_{t}^{2} j_{\parallel} -j_{\parallel} \theta_{t}^{2} -\frac{C \nu \psi^{2}}{2} \\
&\quad -\frac{C \nu \theta^{2}}{2} +C \nu \psi v_{s} +C \nu \theta w_{s}
\end{aligned}
\]


\[
\begin{aligned}
\mathcal E_{c}^{[3]} &=-\Phi_{t}^{2} c j_{\parallel} -c j_{\parallel} \theta_{t}^{2} +\Phi_{t} j_{\parallel} \psi_{t} \theta \\
&\quad -\Phi j_{\parallel} \psi_{t} \theta_{t} -\frac{C \nu \psi^{2} u_{s}}{2} -\frac{C \nu \theta^{2} u_{s}}{2} \\
&\quad +\frac{C \Phi \nu \theta v_{s}}{2} -\frac{C \Phi \nu \psi w_{s}}{2}
\end{aligned}
\]


## Выполненные символические сверки

Все сверки проведены с произвольными символическими коэффициентами и
точными рациональными множителями (1/2, 1/6, 1/24 и т. п.).
Численные значения материала и геометрии не подбирались.

1. Вариация L_(<=4) сопоставлена с независимым разложением неусечённых
   векторных балансов. R, R^T R_s и R^T R_t во втором пути получены
   непосредственно из ряда матричной экспоненты. Материальный остаток
   момента переведён в координатный остаток через P J_r(a)^T,
   P=diag(1,-1,1). Все семь разностей — нулевые многочлены.
2. Развёрнутые потоки сил, матрицы вращательной инерции и кривизн,
   упругий координатный момент и уравнение c совпали с прямыми
   производными L_(<=4).
3. Линейный предел точно возвращает выбранные локальные блоки M-H +
   Тимошенко и Тимошенко + обобщённое кручение.
4. При v=Phi=psi=0 и нулевых производных этих полей три соответствующих
   уравнения обращаются в ноль. Остальные четыре совпадают с независимым
   разложением неусечённой плоской подсистемы.
5. Отражение (v,Phi,psi)->(-v,-Phi,-psi) сохраняет действие. Остатки
   уравнений имеют нужную чётность.
6. Матрица при обобщённых скоростях симметрична.
7. Потоки на границе совпадают с производными V_(<=4) по пространственным
   производным полей. Вращательный поток сопряжён вариации координат
   (Phi,psi,theta), а не непосредственно материальному виртуальному повороту.
8. Локальное тождество энергии проверено как нулевой многочлен:
   d_t(T_(<=4)+V_(<=4))-d_s(Power)-sum_j Q_(j,t) E_(Q_j)=0.
9. Для однородного по s жёсткого поворота с U_s=R(a)e1-e1 и c=0 упругие
   силы и координатные моменты равны нулю до третьей степени.
10. Контроль намеренного пропуска J_r^T не прошёл: все три вращательных
    уравнения получили ненулевые расхождения уже во второй степени.

Это проверки алгебры одной принятой редуцированной модели. Они не являются
сравнением с экспериментом, полной 3D-упругостью, FEM или вычислением
порога выхода колебаний из плоскости. Расчёты спектра, интегрирование
движения и тесты репозитория в этой работе не выполнялись.
