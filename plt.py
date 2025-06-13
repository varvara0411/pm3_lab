## \file plt.py
## \brief Performance comparison visualization for random number generators

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('times.csv')

plt.figure()
plt.plot(df['Size'], df['LCG'], label='LCG')
plt.plot(df['Size'], df['Xorshift32'], label='Xorshift32')
plt.plot(df['Size'], df['MiddleSquare'], label='MiddleSquare')
plt.plot(df['Size'], df['MT19937'], label='std::mt19937')
plt.xlabel('Размер выборки (N)')
plt.ylabel('Время генерации (микросекунды)')

plt.title('Зависимость времени генерации от размера выборки')
plt.legend()
plt.tight_layout()

plt.savefig("time_comparison.png", dpi=300, bbox_inches='tight')
plt.show()

# --- Graph for Xorshift32 and MiddleSquare ---

plt.figure()

plt.plot(df['Size'], df['Xorshift32'], label='Xorshift32')
plt.plot(df['Size'], df['MiddleSquare'], label='MiddleSquare')

plt.xlabel('Размер выборки (N)')
plt.ylabel('Время генерации (микросекунды)')

plt.title('Зависимость времени генерации от размера выборки')
plt.legend()
plt.tight_layout()

plt.savefig("time_comparison_two.png", dpi=300, bbox_inches='tight')
plt.show()