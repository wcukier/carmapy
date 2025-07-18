{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   {% set method_list = methods | reject("equalto", "__init__") | list %}
   {% if method_list %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :signatures: short
   {% for item in methods if item != '__init__' %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes|sort %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% if methods %}
{% for item in methods if item != '__init__' %}
.. automethod:: {{ name }}.{{ item }}
{%- endfor %}

{% endif %}