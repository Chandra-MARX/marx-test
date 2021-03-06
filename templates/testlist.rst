.. _sect-{{ module.__name__ }}:

{{ module.title|default(module.__name__) }}
=============================================================================
{{ module.__doc__ }}

{% for testinstance in testinstances %}

.. _sect-test-{{ testinstance.name }}:

{{ testinstance.title|default(testinstance.name) }}
------------------------------------------------------------------------------

{% if testinstance.obsid %}
:data: `ObsID {{ testinstance.obsid }} <http://cda.harvard.edu/chaser/startViewer.do?menuItem=details&obsid={{testinstance.obsid }}>`_
{% endif %}
:code: :ref:`test-code-{{testinstance.name}}`

{{ testinstance.doc }}

{% for figname, figure in testinstance.figures.items() %}

.. figure:: {{ figpath }}/{{testinstance.name}}_{{ figname }}.*
   :align: center
   :alt: {{ figure.alternative }}
   :scale: {{ figure.scale|default('75%') }}
{% if figure.width %}
   :width: {{ figure.width }}
{% endif %}
{% if figure.height %}
   :height:  {{ figure.height }}
{% endif %}

   {{ figure.caption }}

{% endfor %}

{{ testinstance.summary }}

{% endfor %}
